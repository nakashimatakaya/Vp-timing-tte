################################################################################
# 03_clone_censor_weight.R
#
# Clone-censor-weight framework with steroid-adjusted IPCW.
# Implements the three CCW steps: clone, censor, weight.
#
# Requires: pt_vp6 (with steroid_pre), ts_tbl_data (imputed time-series)
# Outputs:  ccw_tbl_steroid (weighted clone-level table for Cox models)
################################################################################

library(tidyverse)
library(survival)
library(splines)
library(readr)

# ==============================================================================
# 0. Parameters
# ==============================================================================

DAY28_HR  <- 28 * 24
EARLY_END <- 6
P_CLIP    <- 1e-6

# ==============================================================================
# 1. Clone: duplicate each patient into both strategy arms
# ==============================================================================

pt_base_ccw <- pt_vp6 %>%
  filter(!is.na(uid)) %>%
  select(uid, src, any_of(c(
    "age", "sex", "bmi", "wt0",
    "map0", "lact0", "sofa0", "crea0",
    "norepi_equiv0", "mech_vent0",
    "tvp_rel", "death_rel", "death_time", "t0_hr",
    "steroid_pre"
  )))

pt_clone <- bind_rows(
  pt_base_ccw %>% mutate(arm = "ultra"),
  pt_base_ccw %>% mutate(arm = "early")
) %>%
  mutate(
    arm = factor(arm, levels = c("early", "ultra")),
    # Censor: strategy deviation times
    censor_time = case_when(
      arm == "ultra" & (is.na(tvp_rel) | tvp_rel > 3) ~ 3,
      arm == "early" & !is.na(tvp_rel) & tvp_rel <= 3 ~ tvp_rel,
      arm == "early" & (is.na(tvp_rel) | tvp_rel > 6) ~ 6,
      TRUE ~ Inf
    ),
    end_time = pmin(
      coalesce(censor_time, Inf),
      coalesce(death_rel, Inf),
      EARLY_END,
      na.rm = TRUE
    )
  )

# ==============================================================================
# 2. Time-varying covariates for the grace window (0-6h)
# ==============================================================================

ts_tbl_data <- readr::read_rds(
  "/Users/nakashimatakaya/Desktop/SID/Data/Analysis_Data/septic_shock_imputed_df_new"
)
ts_tbl_data <- ts_tbl_data %>%
  mutate(uid = str_c(src, all_id, sep = "::"))

ts_grace_cov <- ts_tbl_data %>%
  inner_join(pt_vp6 %>% select(uid, t0_hr) %>% distinct(), by = "uid") %>%
  mutate(t_rel = time_hr - t0_hr, t = as.integer(round(t_rel))) %>%
  filter(t >= 0, t <= (EARLY_END - 1)) %>%
  group_by(uid, t) %>%
  summarise(
    map = first(map), lact = first(lact),
    norepi_equiv = first(norepi_equiv), sofa = first(sofa),
    mech_vent = first(mech_vent), crea = first(crea),
    epi_rate = first(epi_rate), dopa_rate = first(dopa_rate),
    dobu_rate = first(dobu_rate), phn_rate = first(phn_rate),
    .groups = "drop"
  ) %>%
  mutate(
    mech_vent = if_else(mech_vent == "invasive", "invasive", "noninvasive"),
    mech_vent = replace_na(mech_vent, "noninvasive"),
    mech_vent_inv = as.integer(mech_vent == "invasive")
  ) %>%
  group_by(uid) %>% arrange(t) %>%
  mutate(
    ne0       = first(norepi_equiv),
    ne_from0  = norepi_equiv - ne0,
    ne_delta1 = norepi_equiv - lag(norepi_equiv)
  ) %>%
  ungroup()

# ==============================================================================
# 3. Build grace-window panel data
# ==============================================================================

ccw_grace <- pt_clone %>%
  tidyr::crossing(t = 0:(EARLY_END - 1)) %>%
  filter(t < end_time) %>%
  left_join(ts_grace_cov, by = c("uid", "t")) %>%
  mutate(
    censor_next = is.finite(censor_time) & ((t + 1) == ceiling(censor_time)),
    death_next  = is.finite(death_time) & ((t + 1) == ceiling(death_time / 1)),
    # Missing-indicator imputation
    age_miss  = as.integer(is.na(age)),   age_imp  = replace_na(as.numeric(age), 0),
    bmi_miss  = as.integer(is.na(bmi)),   bmi_imp  = replace_na(as.numeric(bmi), 0),
    wt0_miss  = as.integer(is.na(wt0)),   wt0_imp  = replace_na(as.numeric(wt0), 0),
    map_miss  = as.integer(is.na(map)),   map_imp  = replace_na(as.numeric(map), 0),
    lact_miss = as.integer(is.na(lact)),  lact_imp = replace_na(as.numeric(lact), 0),
    crea_miss = as.integer(is.na(crea)),  crea_imp = replace_na(as.numeric(crea), 0),
    ne_miss   = as.integer(is.na(norepi_equiv)), ne_imp = replace_na(as.numeric(norepi_equiv), 0),
    sofa_miss = as.integer(is.na(sofa)),  sofa_imp = replace_na(as.numeric(sofa), 0),
    ne_from0_miss  = as.integer(is.na(ne_from0)),  ne_from0_imp  = replace_na(as.numeric(ne_from0), 0),
    ne_delta1_miss = as.integer(is.na(ne_delta1)), ne_delta1_imp = replace_na(as.numeric(ne_delta1), 0),
    epi_miss  = as.integer(is.na(epi_rate)),  epi_imp  = replace_na(as.numeric(epi_rate), 0),
    dopa_miss = as.integer(is.na(dopa_rate)), dopa_imp = replace_na(as.numeric(dopa_rate), 0),
    dobu_miss = as.integer(is.na(dobu_rate)), dobu_imp = replace_na(as.numeric(dobu_rate), 0),
    phn_miss  = as.integer(is.na(phn_rate)),  phn_imp  = replace_na(as.numeric(phn_rate), 0),
    steroid_pre = replace_na(as.integer(steroid_pre), 0L)
  )

# ==============================================================================
# 4. Weight: stabilized IPCW with steroid_pre in denominator
# ==============================================================================

fit_ipcw_one_arm <- function(dat) {
  if (sum(dat$censor_next, na.rm = TRUE) < 2) {
    dat$w_ipcw <- 1; dat$w_final <- 1; return(dat)
  }
  sex_term <- if (n_distinct(dat$sex, na.rm = TRUE) >= 2) " + sex" else ""
  src_term <- if (n_distinct(dat$src, na.rm = TRUE) >= 2) " + src" else ""

  f_num <- as.formula(paste0(
    "censor_next ~ ns(t, 3)", sex_term, src_term,
    " + age_imp + age_miss + bmi_imp + bmi_miss + wt0_imp + wt0_miss"))

  f_den <- as.formula(paste0(
    "censor_next ~ ns(t, 3)", sex_term, src_term,
    " + age_imp + age_miss + bmi_imp + bmi_miss + wt0_imp + wt0_miss",
    " + map_imp + map_miss + lact_imp + lact_miss + crea_imp + crea_miss",
    " + ne_imp + ne_miss + ne_from0_imp + ne_from0_miss + ne_delta1_imp + ne_delta1_miss",
    " + epi_imp + epi_miss + dopa_imp + dopa_miss + dobu_imp + dobu_miss + phn_imp + phn_miss",
    " + sofa_imp + sofa_miss + mech_vent_inv",
    " + steroid_pre"))

  fit_num <- glm(f_num, data = dat, family = binomial(), control = glm.control(maxit = 50))
  fit_den <- glm(f_den, data = dat, family = binomial(), control = glm.control(maxit = 50))
  p_num <- pmin(pmax(predict(fit_num, type = "response"), P_CLIP), 1 - P_CLIP)
  p_den <- pmin(pmax(predict(fit_den, type = "response"), P_CLIP), 1 - P_CLIP)

  dat %>%
    mutate(w_ratio = (1 - p_num) / (1 - p_den)) %>%
    group_by(uid, arm) %>% arrange(t) %>%
    mutate(w_cum = cumprod(w_ratio),
           w_ipcw = lag(w_cum, default = 1),
           w_final = last(w_cum)) %>%
    ungroup()
}

cat("Fitting steroid-adjusted IPCW...\n")
w_u <- fit_ipcw_one_arm(ccw_grace %>% filter(arm == "ultra", death_next == FALSE))
w_e <- fit_ipcw_one_arm(ccw_grace %>% filter(arm == "early", death_next == FALSE))
w_grace <- bind_rows(w_u, w_e)

cat("Weight distribution:\n")
print(summary(w_grace$w_final))

# ==============================================================================
# 5. Build final CCW table for outcome analysis
# ==============================================================================

ccw_tbl_steroid <- w_grace %>%
  group_by(uid, arm) %>% arrange(desc(t)) %>% slice(1) %>% ungroup() %>%
  mutate(
    weight = w_final,
    event  = is.finite(death_time) & (death_time <= DAY28_HR),
    tstart = 0,
    tstop  = pmin(coalesce(death_time, Inf), DAY28_HR, na.rm = TRUE)
  )

cap99_s <- as.numeric(quantile(ccw_tbl_steroid$weight, 0.99, na.rm = TRUE))
cap95_s <- as.numeric(quantile(ccw_tbl_steroid$weight, 0.95, na.rm = TRUE))

cat("CCW table: ", nrow(ccw_tbl_steroid), " rows\n")
cat("Cap99:", round(cap99_s, 3), " Cap95:", round(cap95_s, 3), "\n")
