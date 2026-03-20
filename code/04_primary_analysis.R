################################################################################
# 04_primary_analysis.R
#
# Primary outcome: 28-day mortality
# - IPCW Cox (99th/95th percentile truncation)
# - Doubly robust (IPCW + outcome regression)
# - Risk difference and RMST via bootstrap
#
# Requires: ccw_tbl_steroid, ccw_with_outcomes, cap99_s, cap95_s
################################################################################

library(tidyverse)
library(survival)

DAY28_HR <- 28 * 24
B_BOOT   <- 500
SEED     <- 20260319

# ==============================================================================
# 1. Helper
# ==============================================================================

extract_hr <- function(fit, term = "armultra") {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  tibble(HR = ci[1, "exp(coef)"], lower = ci[1, "lower .95"], upper = ci[1, "upper .95"])
}

# ==============================================================================
# 2. Primary Cox models (steroid-adjusted)
# ==============================================================================

run_cox_s <- function(cap_val, label) {
  dat <- ccw_tbl_steroid %>%
    mutate(w = pmax(1e-6, pmin(weight, cap_val)),
           src = factor(src),
           arm = factor(arm, levels = c("early", "ultra")))
  fit <- coxph(Surv(tstart, tstop, event) ~ arm + strata(src),
               data = dat, weights = w, robust = TRUE, cluster = uid)
  extract_hr(fit) %>% mutate(model = label)
}

res_99 <- run_cox_s(cap99_s, "IPCW Cox (99th percentile truncation)")
res_95 <- run_cox_s(cap95_s, "IPCW Cox (95th percentile truncation)")

# Doubly robust
dat_dr <- ccw_with_outcomes %>%
  mutate(w = pmax(1e-6, pmin(weight, cap99_s)),
         arm = factor(arm, levels = c("early", "ultra")),
         src = factor(src))

if (!"mech_vent_inv" %in% names(dat_dr) && "mech_vent0" %in% names(dat_dr)) {
  dat_dr <- dat_dr %>%
    mutate(mech_vent_inv = as.integer(as.character(mech_vent0) == "invasive"))
}

dr_covs <- intersect(
  c("age", "bmi", "wt0", "map0", "lact0", "crea0", "norepi_equiv0",
    "sofa0", "mech_vent_inv", "steroid_pre"),
  names(dat_dr))
for (v in dr_covs) dat_dr[[v]] <- replace_na(as.numeric(dat_dr[[v]]), 0)

dr_formula <- as.formula(paste0(
  "Surv(tstart, tstop, event) ~ arm + strata(src) + ",
  paste(dr_covs, collapse = " + ")))

fit_dr <- tryCatch(
  coxph(dr_formula, data = dat_dr, weights = w, robust = TRUE, cluster = uid),
  error = function(e) {
    coxph(Surv(tstart, tstop, event) ~ arm + strata(src),
          data = dat_dr, weights = w, robust = TRUE, cluster = uid)
  })

res_dr <- extract_hr(fit_dr) %>% mutate(model = "IPCW Cox + outcome regression")
res_primary <- bind_rows(res_99, res_95, res_dr)

cat("=== Primary results (steroid-adjusted) ===\n")
print(res_primary)

# ==============================================================================
# 3. RMST and risk difference via bootstrap
# ==============================================================================

tmp_surv <- ccw_with_outcomes %>%
  mutate(w = pmax(1e-6, pmin(weight, cap99_s)),
         arm = factor(arm, levels = c("early", "ultra")))

sf <- survfit(Surv(tstart, tstop, event) ~ arm, data = tmp_surv, weights = w, id = uid)

calc_rmst_from_sf <- function(sf_obj, strata_name, tau) {
  idx <- which(names(sf_obj$strata) == strata_name)
  if (length(idx) == 0) return(NA_real_)
  s_idx <- if (idx == 1) 1 else (cumsum(sf_obj$strata)[idx - 1] + 1)
  e_idx <- cumsum(sf_obj$strata)[idx]
  t_vec <- sf_obj$time[s_idx:e_idx]
  s_vec <- sf_obj$surv[s_idx:e_idx]
  t_use <- c(0, t_vec[t_vec <= tau], tau)
  s_use <- c(1, s_vec[t_vec <= tau])
  if (length(s_use) < length(t_use)) s_use <- c(s_use, tail(s_use, 1))
  sum(diff(t_use) * head(s_use, -1))
}

sf_df <- tibble(
  time = sf$time, surv = sf$surv,
  arm = rep(names(sf$strata), sf$strata)
) %>% mutate(arm = str_replace(arm, "arm=", ""))

risk_ultra <- 1 - min(sf_df$surv[sf_df$arm == "ultra" & sf_df$time <= DAY28_HR])
risk_early <- 1 - min(sf_df$surv[sf_df$arm == "early" & sf_df$time <= DAY28_HR])
rd_point   <- (risk_ultra - risk_early) * 100

rmst_u_pt    <- calc_rmst_from_sf(sf, "arm=ultra", DAY28_HR)
rmst_e_pt    <- calc_rmst_from_sf(sf, "arm=early", DAY28_HR)
rmst_diff_pt <- (rmst_u_pt - rmst_e_pt) / 24

set.seed(SEED)
boot_fn <- function(data, tau) {
  uids <- unique(data$uid)
  boot_ids <- sample(uids, length(uids), replace = TRUE)
  bd <- map_dfr(seq_along(boot_ids), ~ data %>% filter(uid == boot_ids[.x]) %>% mutate(uid_b = .x))
  sf_b <- tryCatch(
    survfit(Surv(tstart, tstop, event) ~ arm, data = bd, weights = w, id = uid_b),
    error = function(e) NULL)
  if (is.null(sf_b)) return(c(rmst_diff = NA, rd = NA))
  rmst_u <- calc_rmst_from_sf(sf_b, "arm=ultra", tau)
  rmst_e <- calc_rmst_from_sf(sf_b, "arm=early", tau)
  sf_b_df <- tibble(time = sf_b$time, surv = sf_b$surv,
                    arm = rep(names(sf_b$strata), sf_b$strata)) %>%
    mutate(arm = str_replace(arm, "arm=", ""))
  r_u <- 1 - min(sf_b_df$surv[sf_b_df$arm == "ultra" & sf_b_df$time <= tau], na.rm = TRUE)
  r_e <- 1 - min(sf_b_df$surv[sf_b_df$arm == "early" & sf_b_df$time <= tau], na.rm = TRUE)
  c(rmst_diff = (rmst_u - rmst_e) / 24, rd = (r_u - r_e) * 100)
}

cat("Running bootstrap (", B_BOOT, "iterations)...\n")
boot_res <- replicate(B_BOOT, boot_fn(tmp_surv, DAY28_HR))
boot_mat <- t(boot_res)
boot_mat <- boot_mat[complete.cases(boot_mat), ]
rmst_ci  <- quantile(boot_mat[, "rmst_diff"], c(0.025, 0.975))
rd_ci    <- quantile(boot_mat[, "rd"], c(0.025, 0.975))

cat(sprintf("Risk ultra: %.1f%%, Risk early: %.1f%%\n", risk_ultra * 100, risk_early * 100))
cat(sprintf("RD: %.1f pp (%.1f to %.1f)\n", rd_point, rd_ci[1], rd_ci[2]))
cat(sprintf("RMST diff: +%.2f days (%.2f to %.2f)\n", rmst_diff_pt, rmst_ci[1], rmst_ci[2]))
