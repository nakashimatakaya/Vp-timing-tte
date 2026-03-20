################################################################################
# 05_secondary_outcomes.R
#
# Secondary outcomes (weighted Cox, death as competing risk)
# and exploratory subgroup analyses.
#
# Requires: ccw_tbl_steroid, ccw_with_outcomes, pt_vp6, cap99_s
################################################################################

library(tidyverse)
library(survival)

DAY28_HR <- 28 * 24

extract_hr <- function(fit, term = "armultra") {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  tibble(HR = ci[1, "exp(coef)"], lower = ci[1, "lower .95"], upper = ci[1, "upper .95"])
}

# ==============================================================================
# 1. Secondary outcomes
# ==============================================================================

sec_outcomes <- c("aki2_28", "aki3_28", "rrt_28", "crrt_28",
                  "antiarr_28", "shockres_28", "netneg_28")
sec_labels <- c(
  "Acute kidney injury (KDIGO stage \u22652)",
  "Acute kidney injury (KDIGO stage \u22653)",
  "Renal replacement therapy initiation",
  "Continuous renal replacement therapy initiation",
  "Medically treated arrhythmia",
  "Shock resolution",
  "Net negative fluid balance"
)

sec_data <- ccw_tbl_steroid %>%
  left_join(pt_vp6 %>% select(uid, any_of(sec_outcomes)) %>% distinct(), by = "uid")

sec_results <- map2_dfr(sec_outcomes, sec_labels, function(ov, ol) {
  if (!ov %in% names(sec_data))
    return(tibble(outcome = ol, HR = NA_real_, LCL = NA_real_, UCL = NA_real_))
  dat <- sec_data %>%
    filter(!is.na(.data[[ov]])) %>%
    mutate(w = pmax(1e-6, pmin(weight, cap99_s)),
           arm = factor(arm, levels = c("early", "ultra")),
           ev = as.integer(as.logical(.data[[ov]])))
  n_ev <- sum(dat$ev, na.rm = TRUE)
  cat(sprintf("  %s: n=%d, events=%d\n", ol, nrow(dat), n_ev))
  if (n_ev < 5 || n_distinct(dat$arm) < 2)
    return(tibble(outcome = ol, HR = NA_real_, LCL = NA_real_, UCL = NA_real_))
  fit <- tryCatch(
    coxph(Surv(tstart, tstop, ev) ~ arm, data = dat, weights = w, robust = TRUE, cluster = uid),
    error = function(e) NULL)
  if (is.null(fit))
    return(tibble(outcome = ol, HR = NA_real_, LCL = NA_real_, UCL = NA_real_))
  ci <- summary(fit)$conf.int["armultra", , drop = FALSE]
  tibble(outcome = ol, HR = ci[1, "exp(coef)"], LCL = ci[1, "lower .95"], UCL = ci[1, "upper .95"])
})

cat("\n=== Secondary outcomes ===\n")
print(sec_results)

# ==============================================================================
# 2. Subgroup analyses
# ==============================================================================

ccw_sub <- ccw_with_outcomes %>%
  mutate(w = pmax(1e-6, pmin(weight, cap99_s)))

sub_needed <- intersect(c("lact0", "norepi_equiv0", "mech_vent0"), names(pt_vp6))
ccw_sub <- ccw_sub %>%
  select(-any_of(sub_needed)) %>%
  left_join(pt_vp6 %>% select(uid, all_of(sub_needed)) %>% distinct(), by = "uid") %>%
  mutate(
    lact_grp = case_when(lact0 < 2 ~ "<2 mmol/L", lact0 < 4 ~ "2\u2013<4 mmol/L", TRUE ~ "\u22654 mmol/L"),
    ne_grp = if_else(norepi_equiv0 < 0.30, "<0.30 \u00b5g/kg/min", "\u22650.30 \u00b5g/kg/min"),
    mv_grp = if_else(as.character(mech_vent0) == "invasive", "Invasive", "Non-invasive")
  )

fit_sub <- function(dat, grp_var, grp_label) {
  dat %>%
    filter(!is.na(.data[[grp_var]])) %>%
    group_by(level = .data[[grp_var]]) %>%
    group_modify(~ {
      d <- .x %>% mutate(arm = factor(arm, levels = c("early", "ultra")))
      if (n_distinct(d$arm) < 2 || sum(d$event) < 3) return(tibble())
      fit <- tryCatch(
        coxph(Surv(tstart, tstop, event) ~ arm, data = d, weights = w, robust = TRUE, cluster = uid),
        error = function(e) NULL)
      if (is.null(fit)) return(tibble())
      ci <- summary(fit)$conf.int["armultra", , drop = FALSE]
      tibble(HR = ci[1, "exp(coef)"], lower = ci[1, "lower .95"], upper = ci[1, "upper .95"])
    }) %>%
    ungroup() %>%
    mutate(subgroup = grp_label)
}

sub_results <- bind_rows(
  fit_sub(ccw_sub, "lact_grp", "Lactate"),
  fit_sub(ccw_sub, "ne_grp", "NE-equivalent dose"),
  fit_sub(ccw_sub, "mv_grp", "Mechanical ventilation")
)

cat("\n=== Subgroup results ===\n")
print(sub_results)

# ==============================================================================
# 3. Database-specific estimates
# ==============================================================================

db_results <- map_dfr(unique(ccw_with_outcomes$src), function(db) {
  dat_db <- ccw_with_outcomes %>%
    filter(src == db) %>%
    mutate(w = pmax(1e-6, pmin(weight, cap99_s)),
           arm = factor(arm, levels = c("early", "ultra")))
  if (n_distinct(dat_db$arm) < 2 || sum(dat_db$event) < 5)
    return(tibble(database = db, HR = NA_real_, lower = NA_real_, upper = NA_real_))
  fit <- tryCatch(
    coxph(Surv(tstart, tstop, event) ~ arm, data = dat_db, weights = w, robust = TRUE, cluster = uid),
    error = function(e) NULL)
  if (is.null(fit))
    return(tibble(database = db, HR = NA_real_, lower = NA_real_, upper = NA_real_))
  ci <- summary(fit)$conf.int["armultra", , drop = FALSE]
  tibble(database = db, HR = ci[1, "exp(coef)"], lower = ci[1, "lower .95"], upper = ci[1, "upper .95"])
})

cat("\n=== Database-specific results ===\n")
print(db_results)
