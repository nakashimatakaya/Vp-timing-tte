suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(splines)
  library(scales)
  library(flextable)
  library(officer)
  library(rvg)
  library(patchwork)
  library(sandwich)
  library(fs)
})

# Editable figure/table generation for the ICM revision.
# Main analysis: norepinephrine-rate time zero.

ROOT_DIR <- normalizePath(getwd())
DATA_DIR <- Sys.getenv("VP_ANALYSIS_DATA_DIR", unset = "path/to/analysis_data")
REV_DIR <- file.path(DATA_DIR, "revICM")
SHINOZAKI_PT <- Sys.getenv(
  "VP_ELIGIBLE_PATIENT_RDS",
  unset = file.path(DATA_DIR, "pt_all_eligible.rds")
)
TS_IMPUTED_RDS <- Sys.getenv(
  "VP_TS_IMPUTED_RDS",
  unset = file.path(DATA_DIR, "septic_shock_imputed_df_new")
)
FLUID_HOURLY_RDS <- Sys.getenv(
  "VP_FLUID_HOURLY_RDS",
  unset = file.path(DATA_DIR, "_cache_vp_secondary", "fluid_hourly.rds")
)

OUT_DIR <- Sys.getenv("VP_OUTPUT_DIR", unset = file.path(ROOT_DIR, "output", "tables_figures"))
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
QA_DIR <- file.path(OUT_DIR, "qa")
dir_create(c(OUT_DIR, FIG_DIR, TAB_DIR, QA_DIR))

DAY28_HR <- 28 * 24
EARLY_END <- 6
P_CLIP <- 1e-6
SEED <- 20260518
B_BOOT <- as.integer(Sys.getenv("B_BOOT", "500"))

COL_BLACK <- "#111111"
COL_GRID <- "#E7E7E7"
COL_RIBBON <- "#D9D9D9"
FONT_FAMILY <- "Times New Roman"

set.seed(SEED)

fmt_n <- function(x) comma(as.integer(round(x)), accuracy = 1)
fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)
fmt_num <- function(x, digits = 2) ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
fmt_hr <- function(hr, lo, hi) sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
fmt_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    p < 0.10 ~ sprintf("%.3f", p),
    TRUE ~ sprintf("%.2f", p)
  )
}
fmt_med_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"), q[2], q[1], q[3])
}
fmt_count_pct <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  sprintf("%s (%.1f%%)", fmt_n(sum(x)), 100 * mean(x))
}
pp <- function(x) sprintf("%.2f", x)

message("Loading analysis inputs...")
ts_out_rev <- readRDS(file.path(REV_DIR, "vp_ts_out_revICM_20260517.rds"))
ts_imputed <- readRDS(TS_IMPUTED_RDS) %>%
  mutate(uid = paste0(src, "::", as.character(all_id))) %>%
  filter(src %in% c("miiv", "eicu"))

pt_eligible <- readRDS(SHINOZAKI_PT)
if (!"uid" %in% names(pt_eligible)) {
  pt_eligible <- pt_eligible %>% mutate(uid = paste0(src, "::", as.character(all_id)))
}
pt_eligible <- pt_eligible %>% filter(src %in% c("miiv", "eicu"))

comorbidity_tbl <- readRDS(file.path(REV_DIR, "comorbidity_tbl.rds"))
safety_tbl <- readRDS(file.path(REV_DIR, "safety_tbl.rds"))
ne_flags_tbl <- readRDS(file.path(REV_DIR, "ne_flags_tbl.rds"))
t0_neonly_tbl <- readRDS(file.path(REV_DIR, "t0_neonly_tbl.rds"))
inotrope_flags_tbl <- readRDS(file.path(REV_DIR, "inotrope_flags_tbl.rds"))
vp_dose_tbl <- readRDS(file.path(REV_DIR, "vp_dose_tbl.rds"))
vp_interrupt_tbl <- readRDS(file.path(REV_DIR, "vp_interrupt_tbl.rds"))
steroid_timing_tbl_old <- readRDS(file.path(REV_DIR, "steroid_timing_tbl.rds"))

ts_out_secondary <- ts_out_rev %>%
  select(uid, any_of(c(
    "aki2_time", "aki3_time", "aki2_at0", "aki3_at0",
    "rrt_time", "crrt_time", "rrt_pre", "crrt_pre",
    "amio_time", "amio_pre", "net_negative_time", "shock_resolve_time"
  )))

pt_base <- pt_eligible %>%
  distinct(uid, .keep_all = TRUE) %>%
  left_join(comorbidity_tbl, by = "uid") %>%
  left_join(safety_tbl, by = "uid") %>%
  left_join(ne_flags_tbl, by = "uid") %>%
  left_join(t0_neonly_tbl, by = "uid") %>%
  left_join(inotrope_flags_tbl, by = "uid") %>%
  left_join(vp_dose_tbl, by = "uid") %>%
  left_join(vp_interrupt_tbl, by = "uid") %>%
  left_join(steroid_timing_tbl_old, by = "uid") %>%
  left_join(ts_out_secondary, by = "uid") %>%
  mutate(
    across(starts_with("cc_"), ~ coalesce(.x, FALSE)),
    across(starts_with("sf_"), ~ coalesce(.x, FALSE)),
    across(any_of(c(
      "ne_at_t0", "ne_first", "ne_dominant_080", "ne_dominant_090",
      "dobu_at_t0", "dobu_window06", "dobu_pre_t0_24h", "dobu_at_vp",
      "new_dobu_3h", "new_dobu_6h", "new_dobu_24h", "new_epi_24h",
      "interrupt_3h", "sustained_interrupt_6h", "permanent_disc_3h",
      "vp_init_lt003", "vp_init_eq003", "vp_init_gt003", "vp_dose_missing",
      "steroid_pre_t0_24h", "steroid_pre_vp", "steroid_post_vp",
      "aki2_at0", "aki3_at0", "rrt_pre", "crrt_pre", "amio_pre"
    )), ~ coalesce(.x, FALSE))
  )

pt_nee <- pt_base %>%
  mutate(
    t0_hr = as.numeric(t0_hr),
    death_rel = as.numeric(death_time),
    death_time = as.numeric(death_time),
    follow_end_rel = pmin(as.numeric(follow_end_rel), DAY28_HR, na.rm = TRUE),
    tvp_rel = as.numeric(tvp_rel),
    death28 = !is.na(death_time) & is.finite(death_time) &
      death_time >= 0 & death_time <= DAY28_HR
  ) %>%
  filter(!is.na(t0_hr), is.na(tvp_rel) | tvp_rel >= 0)

derive_baseline_at_t0 <- function(pt, ts_imp, t0_col = "t0_hr") {
  ne_rate_col <- intersect(c("norepi_rate", "ne_rate", "norepinephrine_rate"), names(ts_imp))[1]
  ts_imp %>%
    inner_join(pt %>% select(uid, t0_value = all_of(t0_col)), by = "uid") %>%
    mutate(.rel = time_hr - t0_value) %>%
    filter(.rel == 0) %>%
    group_by(uid) %>%
    summarise(
      map0 = first(map),
      lact0 = first(lact),
      crea0 = first(crea),
      sofa0 = first(sofa),
      nee_0 = first(norepi_equiv),
      ne_rate_0 = if (!is.na(ne_rate_col)) first(.data[[ne_rate_col]]) else NA_real_,
      mech_vent0 = first(mech_vent),
      mech_vent0_inv = as.integer(coalesce(first(mech_vent) == "invasive", FALSE)),
      .groups = "drop"
    )
}

derive_steroid_at_t0 <- function(pt, ts_imp, t0_col = "t0_hr") {
  if (!"cort" %in% names(ts_imp)) {
    return(tibble(uid = pt$uid, steroid_pre = as.integer(coalesce(pt$steroid_pre_t0_24h, FALSE))))
  }
  ts_imp %>%
    inner_join(pt %>% select(uid, t0_value = all_of(t0_col)), by = "uid") %>%
    mutate(.rel = time_hr - t0_value, cort_on = coalesce(cort, 0) > 0) %>%
    group_by(uid) %>%
    summarise(steroid_pre = as.integer(any(cort_on & .rel < 0 & .rel >= -24, na.rm = TRUE)),
              .groups = "drop")
}

derive_steroid_fallback <- function(pt) {
  steroid_cols <- intersect(
    c("steroid_pre", "steroid_at_t0", "steroid_pre24h_any", "steroid_pre_t0_24h"),
    names(pt)
  )
  if (!length(steroid_cols)) {
    return(tibble(uid = pt$uid, steroid_pre = 0L))
  }
  pt %>%
    transmute(
      uid,
      steroid_pre = as.integer(if_any(all_of(steroid_cols), ~ coalesce(as.logical(.x), FALSE)))
    )
}

derive_fluid_at_ne_t0 <- function(pt) {
  if (!file.exists(FLUID_HOURLY_RDS)) return(tibble(uid = character()))
  fluid_hourly <- readRDS(FLUID_HOURLY_RDS)
  if (!"uid" %in% names(fluid_hourly)) {
    fluid_hourly <- fluid_hourly %>% mutate(uid = paste0(src, "::", as.character(all_id)))
  }
  if (!"time_hr" %in% names(fluid_hourly) && "hr" %in% names(fluid_hourly)) {
    fluid_hourly <- fluid_hourly %>% rename(time_hr = hr)
  }
  if (!"fluid_balance_ml" %in% names(fluid_hourly) &&
      all(c("fluid_in_ml", "fluid_out_ml") %in% names(fluid_hourly))) {
    fluid_hourly <- fluid_hourly %>%
      mutate(fluid_balance_ml = fluid_in_ml - fluid_out_ml)
  }

  fluid_t0 <- fluid_hourly %>%
    inner_join(pt %>% select(uid, t0_hr), by = "uid") %>%
    mutate(t_rel = time_hr - t0_hr) %>%
    group_by(uid) %>%
    summarise(
      fluid_in_m24_t0 = sum(fluid_in_ml[t_rel >= -24 & t_rel < 0], na.rm = TRUE),
      fluid_in_m6_t0 = sum(fluid_in_ml[t_rel >= -6 & t_rel < 0], na.rm = TRUE),
      fluid_bal_m24_t0 = sum(fluid_balance_ml[t_rel >= -24 & t_rel < 0], na.rm = TRUE),
      fluid_bal_m6_t0 = sum(fluid_balance_ml[t_rel >= -6 & t_rel < 0], na.rm = TRUE),
      .groups = "drop"
    )

  fluid_vp <- fluid_hourly %>%
    inner_join(pt %>% filter(!is.na(tvp_rel)) %>% select(uid, t0_hr, tvp_rel), by = "uid") %>%
    mutate(vp_abs_hr = t0_hr + tvp_rel, t_rel_vp = time_hr - vp_abs_hr) %>%
    group_by(uid) %>%
    summarise(
      fluid_in_m24_vp = sum(fluid_in_ml[t_rel_vp >= -24 & t_rel_vp < 0], na.rm = TRUE),
      fluid_in_m6_vp = sum(fluid_in_ml[t_rel_vp >= -6 & t_rel_vp < 0], na.rm = TRUE),
      fluid_in_p6_vp = sum(fluid_in_ml[t_rel_vp >= 0 & t_rel_vp <= 6], na.rm = TRUE),
      fluid_bal_m24_vp = sum(fluid_balance_ml[t_rel_vp >= -24 & t_rel_vp < 0], na.rm = TRUE),
      fluid_bal_m6_vp = sum(fluid_balance_ml[t_rel_vp >= -6 & t_rel_vp < 0], na.rm = TRUE),
      fluid_bal_p6_vp = sum(fluid_balance_ml[t_rel_vp >= 0 & t_rel_vp <= 6], na.rm = TRUE),
      .groups = "drop"
    )
  fluid_t0 %>% left_join(fluid_vp, by = "uid")
}

shift_time <- function(x, shift) {
  x <- as.numeric(x)
  ifelse(is.na(x) | !is.finite(x), x, x + shift)
}

message("Building norepinephrine-main cohort...")
pt_ne <- pt_base %>%
  filter(!is.na(t0_hr), !is.na(t0_hr_neonly)) %>%
  mutate(
    t0_hr_nee = as.numeric(t0_hr),
    t0_hr = as.numeric(t0_hr_neonly),
    t0_shift = t0_hr_nee - t0_hr,
    tvp_rel = shift_time(tvp_rel, t0_shift),
    death_rel = shift_time(death_time, t0_shift),
    death_time = death_rel,
    follow_end_rel = pmin(shift_time(follow_end_rel, t0_shift), DAY28_HR, na.rm = TRUE),
    across(any_of(c("aki2_time", "aki3_time", "rrt_time", "crrt_time",
                    "amio_time", "net_negative_time", "shock_resolve_time")),
           ~ shift_time(.x, t0_shift))
  ) %>%
  filter(is.na(tvp_rel) | tvp_rel >= 0) %>%
  filter(is.na(death_rel) | !is.finite(death_rel) | death_rel >= 0) %>%
  filter(is.na(follow_end_rel) | follow_end_rel > 0)

baseline_ne <- derive_baseline_at_t0(pt_ne, ts_imputed, "t0_hr")
steroid_ne_raw <- derive_steroid_at_t0(pt_ne, ts_imputed, "t0_hr")
steroid_ne_fallback <- derive_steroid_fallback(pt_ne)
steroid_ne <- if (sum(steroid_ne_raw$steroid_pre == 1L, na.rm = TRUE) > 0) {
  steroid_ne_raw
} else {
  message("Using original steroid flags because the hourly corticosteroid table had no pre-time-zero positives.")
  steroid_ne_fallback
}
fluid_ne <- derive_fluid_at_ne_t0(pt_ne)

pt_ne <- pt_ne %>%
  select(-any_of(c(
    "map0", "lact0", "crea0", "sofa0", "nee_0", "ne_rate_0",
    "mech_vent0", "mech_vent0_inv", "steroid_pre",
    "fluid_in_m24_t0", "fluid_in_m6_t0", "fluid_bal_m24_t0", "fluid_bal_m6_t0",
    "fluid_in_m24_vp", "fluid_in_m6_vp", "fluid_in_p6_vp",
    "fluid_bal_m24_vp", "fluid_bal_m6_vp", "fluid_bal_p6_vp"
  ))) %>%
  left_join(baseline_ne, by = "uid") %>%
  left_join(steroid_ne, by = "uid") %>%
  left_join(fluid_ne, by = "uid") %>%
  mutate(
    death28 = !is.na(death_time) & is.finite(death_time) &
      death_time >= 0 & death_time <= DAY28_HR,
    obs_vp_group = case_when(
      is.na(tvp_rel) ~ "Never VP",
      tvp_rel <= 3 ~ "0-3 h",
      tvp_rel > 3 & tvp_rel <= 6 ~ ">3-6 h",
      tvp_rel > 6 ~ ">6 h",
      TRUE ~ NA_character_
    ),
    obs_vp_group = factor(obs_vp_group, levels = c("0-3 h", ">3-6 h", ">6 h", "Never VP")),
    lact_stratum = case_when(
      is.na(lact0) ~ NA_character_,
      lact0 < 2 ~ "<2 mmol/L",
      lact0 < 4 ~ "2-<4 mmol/L",
      TRUE ~ ">=4 mmol/L"
    ),
    ne_high_035 = case_when(
      is.na(ne_rate_0) ~ NA_integer_,
      ne_rate_0 >= 0.35 ~ 1L,
      TRUE ~ 0L
    ),
    sofa_ge14 = case_when(is.na(sofa0) ~ NA_integer_, sofa0 >= 14 ~ 1L, TRUE ~ 0L),
    aki2_at0 = !is.na(aki2_time) & aki2_time <= 0,
    aki3_at0 = !is.na(aki3_time) & aki3_time <= 0,
    rrt_pre = !is.na(rrt_time) & rrt_time <= 0,
    crrt_pre = !is.na(crrt_time) & crrt_time <= 0,
    amio_pre = !is.na(amio_time) & amio_time <= 0,
    net_negative_at0 = !is.na(net_negative_time) & net_negative_time <= 0
  )

write_csv(pt_ne %>% select(uid, src, t0_hr, t0_hr_nee, t0_shift, tvp_rel, death_rel,
                           follow_end_rel, death28, obs_vp_group),
          file.path(OUT_DIR, "diagnostic_ne_main_cohort.csv"))

safe_min_num <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

make_t0_nad_no_concomitant <- function(df, threshold = 0.25, sustain_k = 3,
                                       sustain_window = 3) {
  cand_tbl <- df %>%
    group_by(uid) %>%
    arrange(time_hr) %>%
    mutate(
      nad_ok = coalesce(norepi_rate, 0) >= threshold,
      other_vasopressor_on = coalesce(epi_rate, 0) > 0 |
        coalesce(dopa_rate, 0) > 0 |
        coalesce(phn_rate, 0) > 0,
      eligible_t0 = nad_ok & !other_vasopressor_on
    ) %>%
    summarise(t0_cand = safe_min_num(time_hr[eligible_t0]), .groups = "drop") %>%
    filter(!is.na(t0_cand))

  df %>%
    select(uid, time_hr, norepi_rate, epi_rate, dopa_rate, phn_rate, adh_rate, vaso_ind) %>%
    inner_join(cand_tbl, by = "uid") %>%
    mutate(
      nad_ok = coalesce(norepi_rate, 0) >= threshold,
      other_vasopressor_on = coalesce(epi_rate, 0) > 0 |
        coalesce(dopa_rate, 0) > 0 |
        coalesce(phn_rate, 0) > 0,
      in_win = time_hr >= t0_cand & time_hr <= (t0_cand + sustain_window),
      ok_in = in_win & nad_ok & !other_vasopressor_on
    ) %>%
    group_by(uid) %>%
    summarise(
      t0_cand = first(t0_cand),
      n_ok = sum(ok_in, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(uid, t0_hr_nad_alone = if_else(n_ok >= sustain_k, t0_cand, NA_real_))
}

build_redefined_t0_cohort <- function(t0_tbl, t0_col) {
  pt_tmp <- pt_base %>%
    left_join(t0_tbl, by = "uid") %>%
    filter(!is.na(t0_hr), !is.na(.data[[t0_col]])) %>%
    mutate(
      t0_hr_nee = as.numeric(t0_hr),
      t0_hr = as.numeric(.data[[t0_col]]),
      t0_shift = t0_hr_nee - t0_hr,
      tvp_rel = shift_time(tvp_rel, t0_shift),
      death_rel = shift_time(death_time, t0_shift),
      death_time = death_rel,
      follow_end_rel = pmin(shift_time(follow_end_rel, t0_shift), DAY28_HR, na.rm = TRUE),
      across(any_of(c("aki2_time", "aki3_time", "rrt_time", "crrt_time",
                      "amio_time", "net_negative_time", "shock_resolve_time")),
             ~ shift_time(.x, t0_shift))
    ) %>%
    filter(is.na(tvp_rel) | tvp_rel >= 0) %>%
    filter(is.na(death_rel) | !is.finite(death_rel) | death_rel >= 0) %>%
    filter(is.na(follow_end_rel) | follow_end_rel > 0)

  baseline_tmp <- derive_baseline_at_t0(pt_tmp, ts_imputed, "t0_hr")
  steroid_tmp_raw <- derive_steroid_at_t0(pt_tmp, ts_imputed, "t0_hr")
  steroid_tmp_fallback <- derive_steroid_fallback(pt_tmp)
  steroid_tmp <- if (sum(steroid_tmp_raw$steroid_pre == 1L, na.rm = TRUE) > 0) {
    steroid_tmp_raw
  } else {
    steroid_tmp_fallback
  }
  fluid_tmp <- derive_fluid_at_ne_t0(pt_tmp)

  pt_tmp %>%
    select(-any_of(c(
      "map0", "lact0", "crea0", "sofa0", "nee_0", "ne_rate_0",
      "mech_vent0", "mech_vent0_inv", "steroid_pre",
      "fluid_in_m24_t0", "fluid_in_m6_t0", "fluid_bal_m24_t0", "fluid_bal_m6_t0",
      "fluid_in_m24_vp", "fluid_in_m6_vp", "fluid_in_p6_vp",
      "fluid_bal_m24_vp", "fluid_bal_m6_vp", "fluid_bal_p6_vp"
    ))) %>%
    left_join(baseline_tmp, by = "uid") %>%
    left_join(steroid_tmp, by = "uid") %>%
    left_join(fluid_tmp, by = "uid") %>%
    mutate(
      death28 = !is.na(death_time) & is.finite(death_time) &
        death_time >= 0 & death_time <= DAY28_HR,
      obs_vp_group = case_when(
        is.na(tvp_rel) ~ "Never VP",
        tvp_rel <= 3 ~ "0-3 h",
        tvp_rel > 3 & tvp_rel <= 6 ~ ">3-6 h",
        tvp_rel > 6 ~ ">6 h",
        TRUE ~ NA_character_
      ),
      obs_vp_group = factor(obs_vp_group, levels = c("0-3 h", ">3-6 h", ">6 h", "Never VP")),
      lact_stratum = case_when(
        is.na(lact0) ~ NA_character_,
        lact0 < 2 ~ "<2 mmol/L",
        lact0 < 4 ~ "2-<4 mmol/L",
        TRUE ~ ">=4 mmol/L"
      ),
      ne_high_035 = case_when(
        is.na(ne_rate_0) ~ NA_integer_,
        ne_rate_0 >= 0.35 ~ 1L,
        TRUE ~ 0L
      ),
      sofa_ge14 = case_when(is.na(sofa0) ~ NA_integer_, sofa0 >= 14 ~ 1L, TRUE ~ 0L),
      aki2_at0 = !is.na(aki2_time) & aki2_time <= 0,
      aki3_at0 = !is.na(aki3_time) & aki3_time <= 0,
      rrt_pre = !is.na(rrt_time) & rrt_time <= 0,
      crrt_pre = !is.na(crrt_time) & crrt_time <= 0,
      amio_pre = !is.na(amio_time) & amio_time <= 0,
      net_negative_at0 = !is.na(net_negative_time) & net_negative_time <= 0
    )
}

message("Building Nad-alone time-zero cohort for Reviewer 1...")
t0_nad_alone_tbl <- make_t0_nad_no_concomitant(ts_imputed)
pt_ne_alone <- build_redefined_t0_cohort(t0_nad_alone_tbl, "t0_hr_nad_alone")
write_csv(pt_ne_alone %>% select(uid, src, t0_hr, t0_hr_nee, t0_shift, tvp_rel,
                                 death_rel, follow_end_rel, death28, obs_vp_group),
          file.path(OUT_DIR, "diagnostic_nad_alone_t0_cohort.csv"))

build_clones <- function(pt_cohort, tvp_var = "tvp_rel", strategy_ultra_h = 3,
                         strategy_early_h = EARLY_END) {
  base <- pt_cohort %>%
    rename(.tvp_rel = all_of(tvp_var)) %>%
    mutate(.tvp_rel = as.numeric(.tvp_rel))

  bind_rows(base %>% mutate(arm = "ultra"),
            base %>% mutate(arm = "early")) %>%
    mutate(
      arm = factor(arm, levels = c("early", "ultra")),
      censor_time = case_when(
        arm == "ultra" & (is.na(.tvp_rel) | .tvp_rel > strategy_ultra_h) ~ strategy_ultra_h,
        arm == "early" & !is.na(.tvp_rel) & .tvp_rel <= strategy_ultra_h ~ .tvp_rel,
        arm == "early" & (is.na(.tvp_rel) | .tvp_rel > strategy_early_h) ~ strategy_early_h,
        TRUE ~ Inf
      ),
      end_time = pmin(censor_time, coalesce(death_rel, Inf), strategy_early_h, na.rm = TRUE),
      tvp_rel = .tvp_rel
    ) %>%
    select(-.tvp_rel)
}

build_grace_panel <- function(pt_clone, ts_imp, t0_var = "t0_hr",
                              early_end_h = EARLY_END) {
  ts_grace_cov <- ts_imp %>%
    inner_join(pt_clone %>% select(uid, all_of(t0_var)) %>% distinct(),
               by = "uid") %>%
    rename(.t0 = all_of(t0_var)) %>%
    mutate(t_rel = time_hr - .t0, t = as.integer(round(t_rel))) %>%
    filter(t >= 0, t <= (early_end_h - 1)) %>%
    group_by(uid, t) %>%
    summarise(
      map = first(map), lact = first(lact), crea = first(crea),
      norepi_equiv = first(norepi_equiv), sofa = first(sofa),
      mech_vent = first(mech_vent),
      epi_rate = first(epi_rate), dopa_rate = first(dopa_rate),
      dobu_rate = first(dobu_rate), phn_rate = first(phn_rate),
      .groups = "drop"
    ) %>%
    mutate(mech_vent_inv = as.integer(coalesce(mech_vent == "invasive", FALSE))) %>%
    group_by(uid) %>%
    arrange(t) %>%
    mutate(
      ne0_panel = first(norepi_equiv),
      ne_from0 = norepi_equiv - ne0_panel,
      ne_delta1 = norepi_equiv - lag(norepi_equiv)
    ) %>%
    ungroup()

  pt_clone %>%
    tidyr::crossing(t = 0:(early_end_h - 1)) %>%
    filter(t < end_time) %>%
    left_join(ts_grace_cov, by = c("uid", "t")) %>%
    mutate(
      censor_next = is.finite(censor_time) & ((t + 1) >= ceiling(censor_time)) &
        (t + 1 <= censor_time + 1e-9),
      death_next = !is.na(death_rel) & is.finite(death_rel) &
        ((t + 1) >= ceiling(death_rel)),
      sex = factor(replace_na(as.character(sex), "Unknown")),
      src = factor(replace_na(as.character(src), "unknown")),
      age_miss = as.integer(is.na(age)), age_imp = replace_na(as.numeric(age), 0),
      bmi_miss = as.integer(is.na(bmi)), bmi_imp = replace_na(as.numeric(bmi), 0),
      wt0_miss = as.integer(is.na(wt0)), wt0_imp = replace_na(as.numeric(wt0), 0),
      map_miss = as.integer(is.na(map)), map_imp = replace_na(as.numeric(map), 0),
      lact_miss = as.integer(is.na(lact)), lact_imp = replace_na(as.numeric(lact), 0),
      crea_miss = as.integer(is.na(crea)), crea_imp = replace_na(as.numeric(crea), 0),
      ne_miss = as.integer(is.na(norepi_equiv)), ne_imp = replace_na(as.numeric(norepi_equiv), 0),
      sofa_miss = as.integer(is.na(sofa)), sofa_imp = replace_na(as.numeric(sofa), 0),
      ne_from0_miss = as.integer(is.na(ne_from0)), ne_from0_imp = replace_na(as.numeric(ne_from0), 0),
      ne_delta1_miss = as.integer(is.na(ne_delta1)), ne_delta1_imp = replace_na(as.numeric(ne_delta1), 0),
      epi_miss = as.integer(is.na(epi_rate)), epi_imp = replace_na(as.numeric(epi_rate), 0),
      dopa_miss = as.integer(is.na(dopa_rate)), dopa_imp = replace_na(as.numeric(dopa_rate), 0),
      dobu_miss = as.integer(is.na(dobu_rate)), dobu_imp = replace_na(as.numeric(dobu_rate), 0),
      phn_miss = as.integer(is.na(phn_rate)), phn_imp = replace_na(as.numeric(phn_rate), 0),
      steroid_pre = replace_na(as.integer(steroid_pre), 0L),
      across(starts_with("cc_"), ~ as.integer(coalesce(.x, FALSE)))
    )
}

fit_ipcw_arm <- function(dat_arm, f_num, f_den) {
  if (sum(dat_arm$censor_next, na.rm = TRUE) < 2) {
    dat_arm$w_ipcw <- 1
    dat_arm$w_cum <- 1
    dat_arm$w_final <- 1
    return(dat_arm)
  }
  clean_prob <- function(p, fallback) {
    p <- as.numeric(p)
    if (!is.finite(fallback)) fallback <- P_CLIP
    p[!is.finite(p)] <- fallback
    pmin(pmax(p, P_CLIP), 1 - P_CLIP)
  }
  fit_num <- glm(f_num, data = dat_arm, family = binomial(), na.action = na.exclude,
                 control = glm.control(maxit = 50))
  fit_den <- glm(f_den, data = dat_arm, family = binomial(), na.action = na.exclude,
                 control = glm.control(maxit = 50))
  fallback_p <- mean(dat_arm$censor_next, na.rm = TRUE)
  p_num <- clean_prob(predict(fit_num, type = "response"), fallback_p)
  p_den <- clean_prob(predict(fit_den, type = "response"), fallback_p)
  if (length(p_num) != nrow(dat_arm) || length(p_den) != nrow(dat_arm)) {
    stop("IPCW prediction length mismatch after covariate preprocessing.")
  }
  dat_arm %>%
    mutate(w_ratio = (1 - p_num) / (1 - p_den)) %>%
    group_by(uid, arm) %>%
    arrange(t) %>%
    mutate(w_cum = cumprod(w_ratio), w_ipcw = lag(w_cum, default = 1), w_final = last(w_cum)) %>%
    ungroup()
}

fit_ipcw_both <- function(ccw_grace, den_terms, num_terms = "ns(t, 3)") {
  f_num <- as.formula(paste0("censor_next ~ ", num_terms))
  f_den <- as.formula(paste0("censor_next ~ ", paste(c(num_terms, den_terms), collapse = " + ")))
  dat <- ccw_grace %>% filter(death_next == FALSE)
  bind_rows(
    fit_ipcw_arm(dat %>% filter(arm == "ultra"), f_num, f_den),
    fit_ipcw_arm(dat %>% filter(arm == "early"), f_num, f_den)
  )
}

build_ccw_terminal <- function(w_grace, day28_hr = DAY28_HR, outcome_time_var = "death_rel") {
  w_grace %>%
    group_by(uid, arm) %>%
    arrange(desc(t)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      out_time_raw = .data[[outcome_time_var]],
      out_time = if_else(is.na(out_time_raw) | is.infinite(out_time_raw), Inf, out_time_raw),
      adm_censor = coalesce(follow_end_rel, day28_hr),
      tstop = pmin(out_time, adm_censor, day28_hr, na.rm = TRUE),
      event = is.finite(out_time_raw) & out_time_raw <= day28_hr,
      weight = w_final,
      tstart = 0
    )
}

run_cox_w <- function(ccw_tbl, cap_val, label, strata_term = "strata(src)", extra_term = NULL) {
  dat <- ccw_tbl %>%
    mutate(w = pmax(P_CLIP, pmin(weight, cap_val)),
           src = factor(src), arm = factor(arm, levels = c("early", "ultra")))
  rhs <- paste(c("arm", strata_term, extra_term), collapse = " + ")
  fit <- coxph(as.formula(paste0("Surv(tstart, tstop, event) ~ ", rhs)),
               data = dat, weights = w, robust = TRUE, cluster = uid)
  ci <- summary(fit)$conf.int["armultra", , drop = FALSE]
  tibble(
    label = label,
    HR = ci[1, "exp(coef)"],
    lower = ci[1, "lower .95"],
    upper = ci[1, "upper .95"],
    n_clone = nrow(dat),
    n_event = sum(dat$event)
  )
}

get_caps <- function(ccw_tbl) {
  list(
    cap99 = as.numeric(quantile(ccw_tbl$weight, 0.99, na.rm = TRUE)),
    cap95 = as.numeric(quantile(ccw_tbl$weight, 0.95, na.rm = TRUE))
  )
}

run_pipeline_full <- function(pt_sub, den_terms, label) {
  pt_clone <- build_clones(pt_sub, "tvp_rel", 3, EARLY_END)
  ccw_grace <- build_grace_panel(pt_clone, ts_imputed, "t0_hr", EARLY_END)
  w_grace <- fit_ipcw_both(ccw_grace, den_terms = den_terms)
  ccw_tbl <- build_ccw_terminal(w_grace, DAY28_HR, "death_rel")
  caps <- get_caps(ccw_tbl)
  list(
    result = run_cox_w(ccw_tbl, caps$cap99, label),
    result95 = run_cox_w(ccw_tbl, caps$cap95, paste0(label, " (95th-percentile truncation)")),
    pt_clone = pt_clone,
    ccw_grace = ccw_grace,
    w_grace = w_grace,
    ccw_tbl = ccw_tbl,
    caps = caps
  )
}

COV_ORIGINAL <- c(
  "sex", "src",
  "age_imp", "age_miss", "bmi_imp", "bmi_miss", "wt0_imp", "wt0_miss",
  "map_imp", "map_miss", "lact_imp", "lact_miss", "crea_imp", "crea_miss",
  "ne_imp", "ne_miss", "ne_from0_imp", "ne_from0_miss", "ne_delta1_imp", "ne_delta1_miss",
  "epi_imp", "epi_miss", "dopa_imp", "dopa_miss", "dobu_imp", "dobu_miss",
  "phn_imp", "phn_miss", "sofa_imp", "sofa_miss", "mech_vent_inv", "steroid_pre"
)
COV_COMORB <- c(
  COV_ORIGINAL,
  "cc_chronic_cardiac", "cc_hypertension", "cc_ckd",
  "cc_diabetes", "cc_malignancy", "cc_hiv", "cc_dementia"
)
COV_NO_INOTROPE <- setdiff(COV_COMORB, c("dobu_imp", "dobu_miss"))

message("Running NE-main CCW/IPCW pipeline...")
main_fit <- run_pipeline_full(pt_ne, COV_COMORB, "Primary: NE time zero")
ccw_tbl_main <- main_fit$ccw_tbl
w_grace_main <- main_fit$w_grace
caps_main <- main_fit$caps

run_database_effect <- function(db_value, db_label) {
  dat <- ccw_tbl_main %>%
    filter(as.character(src) == db_value) %>%
    mutate(
      w = pmax(P_CLIP, pmin(weight, caps_main$cap99)),
      arm = factor(arm, levels = c("early", "ultra"))
    )
  fit <- coxph(Surv(tstart, tstop, event) ~ arm,
               data = dat, weights = w, robust = TRUE, cluster = uid)
  ci <- summary(fit)$conf.int["armultra", , drop = FALSE]
  tibble(
    label = db_label,
    HR = ci[1, "exp(coef)"],
    lower = ci[1, "lower .95"],
    upper = ci[1, "upper .95"],
    n_clone = nrow(dat),
    n_event = sum(dat$event)
  )
}

database_effects <- bind_rows(
  run_database_effect("miiv", "MIMIC-IV"),
  run_database_effect("eicu", "eICU-CRD")
) %>%
  mutate(`HR (95% CI)` = fmt_hr(HR, lower, upper))
write_csv(database_effects, file.path(OUT_DIR, "results_database_specific_ne_main.csv"))

COV_OUTCOME_REG <- intersect(
  c(
    "sex", "age_imp", "age_miss", "bmi_imp", "bmi_miss", "wt0_imp", "wt0_miss",
    "map0", "lact0", "crea0", "sofa0", "ne_rate_0", "nee_0",
    "mech_vent0_inv", "steroid_pre",
    "cc_chronic_cardiac", "cc_hypertension", "cc_ckd",
    "cc_diabetes", "cc_malignancy", "cc_hiv", "cc_dementia"
  ),
  names(ccw_tbl_main)
)
primary_outcome_reg <- run_cox_w(
  ccw_tbl_main,
  caps_main$cap99,
  "IPCW Cox + outcome regression",
  extra_term = paste(COV_OUTCOME_REG, collapse = " + ")
)

message("Running key sensitivity analyses...")
nee_fit <- run_pipeline_full(pt_nee, COV_COMORB, "Sensitivity: NEE time zero")
ne_first_fit <- run_pipeline_full(pt_ne %>% filter(coalesce(ne_first, FALSE)), COV_COMORB,
                                  "Sensitivity: NE-first cohort")
ne_dominant_fit <- run_pipeline_full(pt_ne %>% filter(coalesce(ne_dominant_080, FALSE)), COV_COMORB,
                                     "Sensitivity: NE-dominant at time zero")
no_dobu_t0_fit <- run_pipeline_full(pt_ne %>% filter(!coalesce(dobu_at_t0, FALSE)), COV_COMORB,
                                    "Sensitivity: exclude dobutamine at time zero")
no_dobu_cov_fit <- run_pipeline_full(pt_ne, COV_NO_INOTROPE,
                                     "Sensitivity: censoring model without dobutamine")

derive_t0_vaso_status <- function(pt_sub) {
  pt_lookup <- pt_sub %>%
    select(uid, src_pt = src, t0_hr) %>%
    distinct()

  ts_imputed %>%
    inner_join(pt_lookup, by = "uid") %>%
    mutate(t_rel = time_hr - t0_hr, src = src_pt) %>%
    filter(t_rel == 0) %>%
    group_by(uid, src) %>%
    summarise(
      norepi_rate_t0 = first(replace_na(norepi_rate, 0)),
      epi_rate_t0 = first(replace_na(epi_rate, 0)),
      dopa_rate_t0 = first(replace_na(dopa_rate, 0)),
      phn_rate_t0 = first(replace_na(phn_rate, 0)),
      vasopressin_rate_t0 = first(replace_na(adh_rate, 0)),
      vasopressin_flag_t0 = first(coalesce(vaso_ind, FALSE)),
      dobu_rate_t0 = first(replace_na(dobu_rate, 0)),
      .groups = "drop"
    ) %>%
    mutate(
      nad_on = norepi_rate_t0 > 0,
      vasopressin_on = vasopressin_rate_t0 > 0 | vasopressin_flag_t0,
      other_vasopressor_on = epi_rate_t0 > 0 | dopa_rate_t0 > 0 |
        phn_rate_t0 > 0,
      dobutamine_on = dobu_rate_t0 > 0,
      t0_vasopressor_category = case_when(
        nad_on & !other_vasopressor_on ~ "Nad alone",
        nad_on & other_vasopressor_on ~ "Nad + other vasopressor(s)",
        !nad_on & other_vasopressor_on ~ "No Nad; other vasopressor(s)",
        TRUE ~ "No Nad/other vasopressor recorded"
      )
    )
}

summarise_t0_vaso <- function(status_df, definition_label) {
  levels_cat <- c(
    "Nad alone",
    "Nad + other vasopressor(s)",
    "No Nad; other vasopressor(s)",
    "No Nad/other vasopressor recorded"
  )
  status_df %>%
    mutate(t0_vasopressor_category = factor(t0_vasopressor_category, levels = levels_cat)) %>%
    count(t0_vasopressor_category, name = "n", .drop = FALSE) %>%
    mutate(
      time_zero_definition = definition_label,
      denominator = nrow(status_df),
      percent = if_else(denominator > 0, 100 * n / denominator, NA_real_)
    ) %>%
    select(time_zero_definition, category = t0_vasopressor_category, n, denominator, percent)
}

summarise_t0_vaso_by_db <- function(status_df, definition_label) {
  levels_cat <- c(
    "Nad alone",
    "Nad + other vasopressor(s)",
    "No Nad; other vasopressor(s)",
    "No Nad/other vasopressor recorded"
  )
  status_df %>%
    mutate(
      src = recode(src, miiv = "MIMIC-IV", eicu = "eICU-CRD"),
      t0_vasopressor_category = factor(t0_vasopressor_category, levels = levels_cat)
    ) %>%
    count(src, t0_vasopressor_category, name = "n", .drop = FALSE) %>%
    group_by(src) %>%
    mutate(
      time_zero_definition = definition_label,
      denominator = sum(n),
      percent = if_else(denominator > 0, 100 * n / denominator, NA_real_)
    ) %>%
    ungroup() %>%
    select(time_zero_definition, database = src, category = t0_vasopressor_category,
           n, denominator, percent)
}

summarise_nad_other_combinations <- function(status_df, definition_label) {
  agent_cols <- c("epinephrine", "dopamine", "phenylephrine")
  denominator <- nrow(status_df)
  denominator_nad_other <- sum(status_df$nad_on & status_df$other_vasopressor_on,
                               na.rm = TRUE)

  combo_tbl <- status_df %>%
    filter(nad_on, other_vasopressor_on) %>%
    mutate(
      epinephrine = epi_rate_t0 > 0,
      dopamine = dopa_rate_t0 > 0,
      phenylephrine = phn_rate_t0 > 0,
      n_other_agents = as.integer(epinephrine) + as.integer(dopamine) +
        as.integer(phenylephrine)
    ) %>%
    rowwise() %>%
    mutate(
      other_vasopressor_combination = paste(agent_cols[c_across(all_of(agent_cols))],
                                            collapse = " + ")
    ) %>%
    ungroup() %>%
    count(other_vasopressor_combination, n_other_agents, name = "n") %>%
    arrange(desc(n), other_vasopressor_combination)

  if (!nrow(combo_tbl)) {
    combo_tbl <- tibble(
      other_vasopressor_combination = "No Nad + other vasopressor at time zero",
      n_other_agents = 0L,
      n = 0L
    )
  }

  combo_tbl %>%
    mutate(
      time_zero_definition = definition_label,
      denominator = denominator,
      denominator_nad_plus_other = denominator_nad_other,
      percent_of_definition = if_else(denominator > 0, 100 * n / denominator, NA_real_),
      percent_of_nad_plus_other = if_else(denominator_nad_plus_other > 0,
                                          100 * n / denominator_nad_plus_other,
                                          NA_real_),
      percent_of_definition_label = sprintf("%.1f%%", percent_of_definition),
      percent_of_nad_plus_other_label = if_else(
        is.na(percent_of_nad_plus_other),
        "NA",
        sprintf("%.1f%%", percent_of_nad_plus_other)
      )
    ) %>%
    select(time_zero_definition, other_vasopressor_combination, n_other_agents,
           n, denominator, denominator_nad_plus_other,
           percent_of_definition, percent_of_nad_plus_other,
           percent_of_definition_label, percent_of_nad_plus_other_label)
}

summarise_nad_other_agents <- function(status_df, definition_label) {
  denominator <- nrow(status_df)
  nad_other_df <- status_df %>% filter(nad_on, other_vasopressor_on)
  denominator_nad_other <- nrow(nad_other_df)

  agent_tbl <- tibble(
    other_vasopressor = c("epinephrine", "phenylephrine",
                          "dopamine", "Two or more other vasopressors"),
    n = c(
      sum(nad_other_df$epi_rate_t0 > 0, na.rm = TRUE),
      sum(nad_other_df$phn_rate_t0 > 0, na.rm = TRUE),
      sum(nad_other_df$dopa_rate_t0 > 0, na.rm = TRUE),
      sum((as.integer(nad_other_df$epi_rate_t0 > 0) +
             as.integer(nad_other_df$phn_rate_t0 > 0) +
             as.integer(nad_other_df$dopa_rate_t0 > 0)) >= 2, na.rm = TRUE)
    )
  )

  agent_tbl %>%
    mutate(
      time_zero_definition = definition_label,
      denominator = denominator,
      denominator_nad_plus_other = denominator_nad_other,
      percent_of_definition = if_else(denominator > 0, 100 * n / denominator, NA_real_),
      percent_of_nad_plus_other = if (denominator_nad_other > 0) {
        100 * n / denominator_nad_other
      } else {
        NA_real_
      },
      percent_of_definition_label = sprintf("%.1f%%", percent_of_definition),
      percent_of_nad_plus_other_label = if_else(
        is.na(percent_of_nad_plus_other),
        "NA",
        sprintf("%.1f%%", percent_of_nad_plus_other)
      )
    ) %>%
    select(time_zero_definition, other_vasopressor, n, denominator,
           denominator_nad_plus_other, percent_of_definition,
           percent_of_nad_plus_other, percent_of_definition_label,
           percent_of_nad_plus_other_label)
}

t0_status_nee <- derive_t0_vaso_status(pt_nee)
t0_status_ne_any <- derive_t0_vaso_status(pt_ne)
t0_status_ne_alone <- derive_t0_vaso_status(pt_ne_alone)

message("Running Reviewer 1 NE-only-at-time-zero sensitivity analysis...")
ne_alone_fit <- run_pipeline_full(
  pt_ne_alone,
  COV_COMORB,
  "Nad time zero; no concomitant vasopressor at time zero"
)

reviewer1_t0_vaso_composition <- bind_rows(
  summarise_t0_vaso(t0_status_nee, "NEE >=0.25 mcg/kg/min"),
  summarise_t0_vaso(t0_status_ne_any, "Nad >=0.25 mcg/kg/min; concomitant vasopressors allowed"),
  summarise_t0_vaso(t0_status_ne_alone,
                    "Nad >=0.25 mcg/kg/min; no concomitant vasopressors at time zero")
) %>%
  mutate(percent_label = sprintf("%.1f%%", percent))
write_csv(reviewer1_t0_vaso_composition,
          file.path(OUT_DIR, "reviewer1_timezero_vasopressor_composition_20260519.csv"))

reviewer1_t0_vaso_composition_by_db <- bind_rows(
  summarise_t0_vaso_by_db(t0_status_nee, "NEE >=0.25 mcg/kg/min"),
  summarise_t0_vaso_by_db(t0_status_ne_any, "Nad >=0.25 mcg/kg/min; concomitant vasopressors allowed"),
  summarise_t0_vaso_by_db(t0_status_ne_alone,
                          "Nad >=0.25 mcg/kg/min; no concomitant vasopressors at time zero")
) %>%
  mutate(percent_label = sprintf("%.1f%%", percent))
write_csv(reviewer1_t0_vaso_composition_by_db,
          file.path(OUT_DIR, "reviewer1_timezero_vasopressor_composition_by_database_20260519.csv"))

reviewer1_nad_other_combination_breakdown <- bind_rows(
  summarise_nad_other_combinations(t0_status_nee, "NEE >=0.25 mcg/kg/min"),
  summarise_nad_other_combinations(
    t0_status_ne_any,
    "Nad >=0.25 mcg/kg/min; concomitant vasopressors allowed"
  ),
  summarise_nad_other_combinations(
    t0_status_ne_alone,
    "Nad >=0.25 mcg/kg/min; no concomitant vasopressors at time zero"
  )
)
write_csv(reviewer1_nad_other_combination_breakdown,
          file.path(OUT_DIR, "reviewer1_nad_other_combination_breakdown_20260519.csv"))

reviewer1_nad_other_agent_breakdown <- bind_rows(
  summarise_nad_other_agents(t0_status_nee, "NEE >=0.25 mcg/kg/min"),
  summarise_nad_other_agents(
    t0_status_ne_any,
    "Nad >=0.25 mcg/kg/min; concomitant vasopressors allowed"
  ),
  summarise_nad_other_agents(
    t0_status_ne_alone,
    "Nad >=0.25 mcg/kg/min; no concomitant vasopressors at time zero"
  )
)
write_csv(reviewer1_nad_other_agent_breakdown,
          file.path(OUT_DIR, "reviewer1_nad_other_agent_breakdown_20260519.csv"))

primary_sens <- bind_rows(
  main_fit$result,
  main_fit$result95,
  primary_outcome_reg,
  nee_fit$result,
  ne_first_fit$result,
  ne_dominant_fit$result,
  no_dobu_t0_fit$result,
  no_dobu_cov_fit$result
) %>%
  mutate(HR_95CI = fmt_hr(HR, lower, upper))
write_csv(primary_sens, file.path(OUT_DIR, "results_primary_sensitivity_ne_main.csv"))

or11_sensitivity_results <- bind_rows(
  main_fit$result %>% mutate(label = "Primary: NE time zero"),
  main_fit$result95 %>% mutate(label = "Primary: 95th percentile truncation"),
  database_effects %>%
    select(label, HR, lower, upper, n_clone, n_event) %>%
    mutate(label = paste0("Database: ", label)),
  nee_fit$result %>% mutate(label = "Sensitivity: NEE time zero"),
  ne_alone_fit$result %>% mutate(label = "Sensitivity: Nad-only time zero"),
  no_dobu_t0_fit$result %>% mutate(label = "Sensitivity: excluding dobutamine at time zero")
) %>%
  mutate(HR_95CI = fmt_hr(HR, lower, upper))
write_csv(or11_sensitivity_results,
          file.path(OUT_DIR, "results_online_resource11_sensitivity_forest.csv"))

weighted_km <- function(dat, cap_val, t_max = DAY28_HR) {
  d <- dat %>%
    mutate(w = pmax(P_CLIP, pmin(weight, cap_val)),
           arm = factor(arm, levels = c("early", "ultra")))
  map_dfr(c("early", "ultra"), function(a) {
    s <- d %>% filter(arm == a)
    fit <- survfit(Surv(tstart, tstop, event) ~ 1, data = s, weights = s$w)
    times <- fit$time[fit$time <= t_max]
    surv <- fit$surv[fit$time <= t_max]
    if (!length(times) || times[1] > 0) {
      times <- c(0, times)
      surv <- c(1, surv)
    }
    if (max(times, na.rm = TRUE) < t_max) {
      times <- c(times, t_max)
      surv <- c(surv, tail(surv, 1))
    }
    tibble(arm = a, time = times, cuminc = 1 - surv)
  })
}

get_riskdiff <- function(km_df, t_max = DAY28_HR) {
  risk_at <- function(a) {
    km_df %>%
      filter(arm == a, time <= t_max) %>%
      arrange(time) %>%
      slice_tail(n = 1) %>%
      pull(cuminc)
  }
  r_u <- risk_at("ultra")
  r_e <- risk_at("early")
  list(r_u = r_u, r_e = r_e, rd = r_u - r_e)
}

get_rmst_diff <- function(km_df, t_max = DAY28_HR) {
  rmst_arm <- function(df) {
    df <- df %>% arrange(time) %>% filter(time <= t_max)
    if (tail(df$time, 1) < t_max) {
      df <- bind_rows(df, tibble(arm = df$arm[1], time = t_max, cuminc = tail(df$cuminc, 1)))
    }
    surv <- 1 - df$cuminc
    sum(diff(df$time) * head(surv, -1)) / 24
  }
  rmst_arm(km_df %>% filter(arm == "ultra")) -
    rmst_arm(km_df %>% filter(arm == "early"))
}

km_main <- weighted_km(ccw_tbl_main, caps_main$cap99)
pt_rd <- get_riskdiff(km_main)
pt_rmst <- get_rmst_diff(km_main)

boot_rd_rmst <- function(ccw_tbl, cap_val, B = B_BOOT) {
  if (B <= 0) return(tibble(b = integer(), rd = numeric(), rmst = numeric()))
  uids <- unique(ccw_tbl$uid)
  map_dfr(seq_len(B), function(b) {
    sampled <- sample(uids, length(uids), replace = TRUE)
    bdat <- tibble(uid = sampled) %>%
      inner_join(ccw_tbl, by = "uid", relationship = "many-to-many")
    km_b <- weighted_km(bdat, cap_val)
    tibble(b = b, rd = get_riskdiff(km_b)$rd, rmst = get_rmst_diff(km_b))
  })
}

message("Bootstrapping risk difference / RMST: B = ", B_BOOT)
boot_out <- boot_rd_rmst(ccw_tbl_main, caps_main$cap99, B_BOOT)
ci_rd <- quantile(boot_out$rd, c(0.025, 0.975), na.rm = TRUE)
ci_rmst <- quantile(boot_out$rmst, c(0.025, 0.975), na.rm = TRUE)
rd_rmst_results <- tibble(
  metric = c("28-d risk (ultra)", "28-d risk (early)", "RD (pp)", "RMST diff (days)"),
  estimate = c(pt_rd$r_u, pt_rd$r_e, pt_rd$rd * 100, pt_rmst),
  lower = c(NA_real_, NA_real_, ci_rd[1] * 100, ci_rmst[1]),
  upper = c(NA_real_, NA_real_, ci_rd[2] * 100, ci_rmst[2])
)
write_csv(rd_rmst_results, file.path(OUT_DIR, "results_primary_rd_rmst_ne_main.csv"))

summarise_reviewer1_effect <- function(definition_label, pt_sub, fit_obj) {
  km <- weighted_km(fit_obj$ccw_tbl, fit_obj$caps$cap99)
  rd <- get_riskdiff(km)
  rmst <- get_rmst_diff(km)
  res <- fit_obj$result %>% slice(1)

  tibble(
    time_zero_definition = definition_label,
    n_patients = nrow(pt_sub),
    n_mimic = sum(pt_sub$src == "miiv", na.rm = TRUE),
    n_eicu = sum(pt_sub$src == "eicu", na.rm = TRUE),
    n_adherent_ultra_0_3h = sum(!is.na(pt_sub$tvp_rel) & pt_sub$tvp_rel <= 3),
    n_adherent_early_gt3_6h = sum(!is.na(pt_sub$tvp_rel) & pt_sub$tvp_rel > 3 &
                                    pt_sub$tvp_rel <= 6),
    HR = res$HR,
    lower = res$lower,
    upper = res$upper,
    HR_95CI = fmt_hr(res$HR, res$lower, res$upper),
    risk_ultra = rd$r_u,
    risk_early = rd$r_e,
    risk_difference_pp = rd$rd * 100,
    rmst_difference_days = rmst
  )
}

reviewer1_timezero_effects <- bind_rows(
  summarise_reviewer1_effect("NEE >=0.25 mcg/kg/min", pt_nee, nee_fit),
  summarise_reviewer1_effect("Nad >=0.25 mcg/kg/min; concomitant vasopressors allowed",
                             pt_ne, main_fit),
  summarise_reviewer1_effect("Nad >=0.25 mcg/kg/min; no concomitant vasopressors at time zero",
                             pt_ne_alone, ne_alone_fit)
) %>%
  mutate(
    risk_ultra_label = sprintf("%.1f%%", 100 * risk_ultra),
    risk_early_label = sprintf("%.1f%%", 100 * risk_early),
    risk_difference_label = sprintf("%.1f pp", risk_difference_pp),
    rmst_difference_label = sprintf("%.2f days", rmst_difference_days)
  )
write_csv(reviewer1_timezero_effects,
          file.path(OUT_DIR, "reviewer1_timezero_definition_effects_20260519.csv"))

reviewer1_timezero_summary <- reviewer1_t0_vaso_composition %>%
  filter(category %in% c("Nad alone", "Nad + other vasopressor(s)",
                         "No Nad; other vasopressor(s)")) %>%
  select(time_zero_definition, category, n, percent_label) %>%
  pivot_wider(
    names_from = category,
    values_from = c(n, percent_label),
    names_sep = "__"
  ) %>%
  left_join(
    reviewer1_timezero_effects %>%
      select(time_zero_definition, n_patients, n_mimic, n_eicu,
             n_adherent_ultra_0_3h, n_adherent_early_gt3_6h,
             HR_95CI, risk_ultra_label, risk_early_label,
             risk_difference_label, rmst_difference_label),
    by = "time_zero_definition"
  )
write_csv(reviewer1_timezero_summary,
          file.path(OUT_DIR, "reviewer1_timezero_definition_summary_20260519.csv"))

run_finegray_one <- function(w_grace_in, time_var, at0_var, label, cap_val = caps_main$cap99) {
  attach <- pt_ne %>%
    transmute(uid,
              out_time_secondary = .data[[time_var]],
              at0_secondary = coalesce(.data[[at0_var]], FALSE))
  base <- w_grace_in %>%
    group_by(uid, arm) %>%
    arrange(desc(t)) %>%
    slice(1) %>%
    ungroup() %>%
    left_join(attach, by = "uid") %>%
    mutate(
      event_time = if_else(!is.na(out_time_secondary) & out_time_secondary > 0 &
                             !at0_secondary, out_time_secondary, Inf),
      death_time_fg = if_else(!is.na(death_rel) & is.finite(death_rel) & death_rel >= 0,
                              death_rel, Inf),
      adm_censor = coalesce(follow_end_rel, DAY28_HR),
      tstop = pmin(event_time, death_time_fg, adm_censor, DAY28_HR, na.rm = TRUE),
      status = case_when(
        is.finite(event_time) & event_time <= tstop & event_time <= DAY28_HR ~ "event",
        is.finite(death_time_fg) & death_time_fg <= tstop & death_time_fg <= DAY28_HR ~ "death",
        TRUE ~ "censor"
      ),
      status_f = factor(status, levels = c("censor", "event", "death")),
      w = pmax(P_CLIP, pmin(w_final, cap_val)),
      arm = factor(arm, levels = c("early", "ultra")),
      src = factor(src)
    ) %>%
    filter(tstop > 0)

  if (sum(base$status == "event") < 2) {
    return(tibble(label = label, HR = NA_real_, lower = NA_real_, upper = NA_real_,
                  n_clone = nrow(base), n_event = sum(base$status == "event")))
  }
  fg <- finegray(Surv(tstop, status_f) ~ ., data = base, etype = "event", weights = w)
  fit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ arm + strata(src),
               data = fg, weights = fgwt, robust = TRUE, cluster = uid)
  ci <- summary(fit)$conf.int["armultra", , drop = FALSE]
  tibble(label = label, HR = ci[1, "exp(coef)"], lower = ci[1, "lower .95"],
         upper = ci[1, "upper .95"], n_clone = nrow(base),
         n_event = sum(base$status == "event"))
}

secondary_finegray <- bind_rows(
  run_finegray_one(w_grace_main, "aki2_time", "aki2_at0", "AKI (KDIGO stage >=2)"),
  run_finegray_one(w_grace_main, "aki3_time", "aki3_at0", "AKI (KDIGO stage >=3)"),
  run_finegray_one(w_grace_main, "rrt_time", "rrt_pre", "RRT initiation"),
  run_finegray_one(w_grace_main, "crrt_time", "crrt_pre", "CRRT initiation"),
  run_finegray_one(w_grace_main, "amio_time", "amio_pre", "Medically treated arrhythmia"),
  run_finegray_one(w_grace_main, "net_negative_time", "net_negative_at0", "Net negative fluid balance")
) %>%
  mutate(HR_95CI = fmt_hr(HR, lower, upper))
write_csv(secondary_finegray, file.path(OUT_DIR, "results_secondary_finegray_ne_main.csv"))

run_subgroup_interaction <- function(ccw_tbl, sub_var, label) {
  dat <- ccw_tbl %>%
    select(-any_of(sub_var)) %>%
    left_join(pt_ne %>% select(uid, all_of(sub_var)), by = "uid") %>%
    mutate(w = pmax(P_CLIP, pmin(weight, caps_main$cap99)),
           arm = factor(arm, levels = c("early", "ultra")),
           src = factor(src), .sub = .data[[sub_var]]) %>%
    filter(!is.na(.sub))
  if (is.logical(dat$.sub)) dat$.sub <- as.integer(dat$.sub)
  if (!is.numeric(dat$.sub)) dat$.sub <- factor(dat$.sub)
  fit <- coxph(Surv(tstart, tstop, event) ~ arm * .sub + strata(src),
               data = dat, weights = w, robust = TRUE, cluster = uid)
  coef_int <- grep(":", names(coef(fit)), value = TRUE)
  vc <- vcov(fit)
  est <- coef(fit)[coef_int]
  V <- vc[coef_int, coef_int, drop = FALSE]
  p_int <- tryCatch({
    if (length(est) == 1) {
      2 * pnorm(-abs(est / sqrt(V)))
    } else {
      chi <- as.numeric(t(est) %*% solve(V) %*% est)
      pchisq(chi, df = length(est), lower.tail = FALSE)
    }
  }, error = function(e) NA_real_)
  tibble(label = label, sub_var = sub_var, interaction_p = as.numeric(p_int))
}

run_within_stratum <- function(ccw_tbl, stratum_var) {
  ccw_tbl %>%
    select(-any_of(stratum_var)) %>%
    left_join(pt_ne %>% select(uid, all_of(stratum_var)), by = "uid") %>%
    filter(!is.na(.data[[stratum_var]])) %>%
    group_by(.stratum = .data[[stratum_var]]) %>%
    group_modify(~ {
      d <- .x %>%
        mutate(w = pmax(P_CLIP, pmin(weight, caps_main$cap99)),
               arm = factor(arm, levels = c("early", "ultra")),
               src = factor(src))
      fit <- if (length(unique(d$src)) < 2) {
        coxph(Surv(tstart, tstop, event) ~ arm, data = d, weights = w,
              robust = TRUE, cluster = uid)
      } else {
        coxph(Surv(tstart, tstop, event) ~ arm + strata(src), data = d,
              weights = w, robust = TRUE, cluster = uid)
      }
      ci <- summary(fit)$conf.int["armultra", , drop = FALSE]
      tibble(HR = ci[1, "exp(coef)"], lower = ci[1, "lower .95"],
             upper = ci[1, "upper .95"], n_clone = nrow(d), n_event = sum(d$event))
    }) %>%
    ungroup() %>%
    mutate(.stratum = as.character(.stratum), stratum_var = stratum_var)
}

subgroup_interactions <- bind_rows(
  run_subgroup_interaction(ccw_tbl_main, "lact_stratum", "Lactate at time zero"),
  run_subgroup_interaction(ccw_tbl_main, "ne_high_035", "Norepinephrine dose >=0.35"),
  run_subgroup_interaction(ccw_tbl_main, "steroid_pre", "Baseline steroid"),
  run_subgroup_interaction(ccw_tbl_main, "mech_vent0_inv", "Invasive mechanical ventilation"),
  run_subgroup_interaction(ccw_tbl_main, "sofa_ge14", "SOFA score at time zero")
)
subgroup_within <- bind_rows(
  run_within_stratum(ccw_tbl_main, "lact_stratum"),
  run_within_stratum(ccw_tbl_main, "ne_high_035"),
  run_within_stratum(ccw_tbl_main, "steroid_pre"),
  run_within_stratum(ccw_tbl_main, "mech_vent0_inv"),
  run_within_stratum(ccw_tbl_main, "sofa_ge14")
) %>%
  mutate(HR_95CI = fmt_hr(HR, lower, upper))
write_csv(subgroup_interactions, file.path(OUT_DIR, "results_subgroup_interaction_ne_main.csv"))
write_csv(subgroup_within, file.path(OUT_DIR, "results_subgroup_within_ne_main.csv"))

theme_icm <- function(base_size = 18) {
  theme_minimal(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      text = element_text(color = COL_BLACK, family = FONT_FAMILY, face = "bold"),
      axis.text = element_text(color = "#4A4A4A", face = "bold"),
      axis.title = element_text(color = COL_BLACK, face = "bold"),
      panel.grid.major = element_line(color = COL_GRID, linewidth = 0.55),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 8),
      plot.margin = margin(16, 24, 16, 16)
    )
}

make_figure1 <- function() {
  n_total <- nrow(pt_ne)
  n_miiv <- sum(pt_ne$src == "miiv")
  n_eicu <- sum(pt_ne$src == "eicu")
  n_ultra <- sum(!is.na(pt_ne$tvp_rel) & pt_ne$tvp_rel <= 3)
  n_early <- sum(!is.na(pt_ne$tvp_rel) & pt_ne$tvp_rel > 3 & pt_ne$tvp_rel <= 6)
  n_cens_ultra <- n_total - n_ultra
  n_cens_early <- n_total - n_early
  ge <- "\u2265"
  mu <- "\u00b5"
  box_df <- function(id, x, y, w, h, label) {
    tibble(id = id, x = x, y = y, w = w, h = h, label = label,
           xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2)
  }
  boxes <- bind_rows(
    box_df("source", 0.50, 0.93, 0.50, 0.09,
           paste0("Adults with septic shock in MIMIC-IV and eICU-CRD\n",
                  "N = 21,587\n(MIMIC-IV 14,531; eICU-CRD 7,056)")),
    box_df("eligible", 0.50, 0.785, 0.42, 0.12,
           paste0("Eligible at time zero\nNE dose ", ge, "0.25 ", mu, "g/kg/min, no prior VP\n",
                  "n = ", fmt_n(n_total), "\n(MIMIC-IV ", fmt_n(n_miiv),
                  "; eICU-CRD ", fmt_n(n_eicu), ")")),
    box_df("cloned", 0.50, 0.63, 0.46, 0.095,
           paste0("Cloned patient records\nEach patient represented once under each strategy\n",
                  "n = ", fmt_n(n_total), " x 2 = ", fmt_n(2 * n_total), " clones")),
    box_df("assign_ultra", 0.28, 0.455, 0.30, 0.085,
           paste0("Assigned to ultra-early strategy\n0-3 h\nn = ", fmt_n(n_total), " clones")),
    box_df("assign_early", 0.72, 0.455, 0.30, 0.085,
           paste0("Assigned to early strategy\n>3-6 h\nn = ", fmt_n(n_total), " clones")),
    box_df("cens_ultra", 0.12, 0.30, 0.22, 0.075,
           paste0("Artificial censoring due\nto nonadherence\nn = ", fmt_n(n_cens_ultra))),
    box_df("cens_early", 0.88, 0.30, 0.22, 0.075,
           paste0("Artificial censoring due\nto nonadherence\nn = ", fmt_n(n_cens_early))),
    box_df("adh_ultra", 0.28, 0.14, 0.30, 0.09,
           paste0("Adhered to assigned strategy\nVP initiated within 0-3 h\nn = ",
                  fmt_n(n_ultra), " (", sprintf("%.1f", 100 * n_ultra / n_total), "%)")),
    box_df("adh_early", 0.72, 0.14, 0.30, 0.09,
           paste0("Adhered to assigned strategy\nVP initiated within >3-6 h\nn = ",
                  fmt_n(n_early), " (", sprintf("%.1f", 100 * n_early / n_total), "%)"))
  )
  connectors <- tribble(
    ~x, ~y, ~xend, ~yend,
    0.50, 0.582, 0.50, 0.535,
    0.50, 0.535, 0.28, 0.535,
    0.50, 0.535, 0.72, 0.535
  )
  arrows <- tribble(
    ~x, ~y, ~xend, ~yend,
    0.50, 0.885, 0.50, 0.850,
    0.50, 0.725, 0.50, 0.680,
    0.28, 0.535, 0.28, 0.505,
    0.72, 0.535, 0.72, 0.505,
    0.28, 0.410, 0.28, 0.190,
    0.72, 0.410, 0.72, 0.190,
    0.28, 0.300, 0.235, 0.300,
    0.72, 0.300, 0.765, 0.300
  )

  ggplot() +
    geom_segment(data = connectors, aes(x, y, xend = xend, yend = yend),
                 linewidth = 0.7, color = COL_BLACK) +
    geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.14, "inches"), type = "closed"),
                 linewidth = 0.7, color = COL_BLACK) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", color = COL_BLACK, linewidth = 0.65) +
    geom_text(data = boxes, aes(x, y, label = label),
              size = 4.05, lineheight = 0.92, family = FONT_FAMILY,
              fontface = "plain", color = COL_BLACK) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

fig1 <- make_figure1()

rd_text <- sprintf(
  "Day-28 mortality risk\n0-3 h: %.1f%% vs. >3-6 h: %.1f%%\n\nRisk difference: %.1f pp (%.1f to %.1f)\n\nRMST difference: +%.2f days (%.2f to %.2f)",
  100 * pt_rd$r_u, 100 * pt_rd$r_e,
  pt_rd$rd * 100, ci_rd[1] * 100, ci_rd[2] * 100,
  pt_rmst, ci_rmst[1], ci_rmst[2]
)

fig2 <- km_main %>%
  mutate(arm = factor(arm, levels = c("ultra", "early"),
                      labels = c("0-3 h (ultra-early)", ">3-6 h (early)"))) %>%
  ggplot(aes(x = time / 24, y = cuminc * 100, linetype = arm)) +
  geom_step(linewidth = 1.25, color = COL_BLACK) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_x_continuous(breaks = seq(0, 28, 7), limits = c(0, 28)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  annotate("label", x = 14.2, y = 7, label = rd_text,
           hjust = 0, vjust = 0, size = 5.1, family = FONT_FAMILY,
           fontface = "bold", fill = "white", color = COL_BLACK,
           linewidth = 0, lineheight = 1.02) +
  labs(x = "Days after time zero", y = "Cumulative incidence of death (%)",
       linetype = NULL) +
  theme_icm(20) +
  theme(
    legend.position = c(0.78, 0.92),
    legend.text = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    axis.title.x = element_text(size = 27, margin = margin(t = 16)),
    axis.title.y = element_text(size = 25, margin = margin(r = 14))
  )

make_forest <- function(df, title = NULL, xlim = c(0.55, 1.18), breaks = NULL) {
  plot_df <- df %>%
    mutate(row = rev(row_number()), text = fmt_hr(HR, lower, upper))
  if (is.null(breaks)) breaks <- pretty(xlim, 5)
  ggplot(plot_df, aes(y = row, x = HR)) +
    geom_vline(xintercept = 1, linetype = "22", linewidth = 0.8, color = "#555555") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                  width = 0.13, linewidth = 0.95, color = COL_BLACK) +
    geom_point(size = 4.6, color = COL_BLACK) +
    geom_text(aes(x = xlim[2] - 0.015, label = text),
              hjust = 1, size = 5.0, family = FONT_FAMILY, fontface = "bold") +
    scale_y_continuous(breaks = plot_df$row, labels = plot_df$label,
                       expand = expansion(add = 0.6)) +
    scale_x_continuous(limits = xlim, breaks = breaks) +
    labs(title = title, x = "Hazard ratio (ultra-early vs early)", y = NULL) +
    theme_icm(18) +
    theme(axis.text.y = element_text(size = 16, face = "bold", hjust = 1),
          axis.title.x = element_text(size = 21, margin = margin(t = 10)),
          plot.title = element_text(size = 25, face = "bold"))
}

fig3_df <- primary_sens %>%
  filter(label %in% c(
    "Primary: NE time zero",
    "Primary: NE time zero (95th-percentile truncation)"
  )) %>%
  mutate(
    label = recode(
      label,
      "Primary: NE time zero" = "IPCW Cox (99th percentile truncation)",
      "Primary: NE time zero (95th-percentile truncation)" = "IPCW Cox (95th percentile truncation)"
    ),
    order = match(
      label,
      c(
        "IPCW Cox (99th percentile truncation)",
        "IPCW Cox (95th percentile truncation)"
      )
    )
  ) %>%
  arrange(order) %>%
  transmute(label, HR, lower, upper, row = rev(row_number()), text = fmt_hr(HR, lower, upper))

fig3 <- ggplot(fig3_df, aes(y = row, x = HR)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.9, color = "#555555") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                width = 0.11, linewidth = 1.05, color = COL_BLACK) +
  geom_point(size = 5.1, color = COL_BLACK) +
  geom_text(aes(label = text), nudge_y = 0.23, size = 5.1,
            family = FONT_FAMILY, fontface = "bold", color = "#333333") +
  scale_y_continuous(breaks = fig3_df$row, labels = fig3_df$label,
                     limits = c(0.35, 2.65), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.45, 1.08), breaks = seq(0.5, 1.0, 0.1)) +
  labs(title = "Primary outcome", x = "Hazard ratio (ultra-early vs early)", y = NULL) +
  theme_icm(22) +
  theme(
    panel.border = element_rect(fill = NA, color = "#333333", linewidth = 1.0),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 16.5, face = "bold", hjust = 1),
    axis.title.x = element_text(size = 20.5, margin = margin(t = 10)),
    plot.title = element_text(size = 25, face = "bold", hjust = 0),
    plot.margin = margin(18, 28, 18, 18)
  )

fig3a_df <- primary_sens %>%
  filter(label != "IPCW Cox + outcome regression") %>%
  transmute(label = recode(label,
                           "Primary: NE time zero" = "Primary: NE time zero",
                           "Primary: NE time zero (95th-percentile truncation)" = "Primary: 95th percentile truncation",
                           "IPCW Cox + outcome regression" = "IPCW Cox + outcome regression",
                           "Sensitivity: NEE time zero" = "Sensitivity: NEE time zero",
                           "Sensitivity: NE-first cohort" = "Sensitivity: NE-first cohort",
                           "Sensitivity: NE-dominant at time zero" = "Sensitivity: NE-dominant at time zero",
                           "Sensitivity: exclude dobutamine at time zero" = "Sensitivity: exclude dobutamine at time zero",
                           "Sensitivity: censoring model without dobutamine" = "Sensitivity: no dobutamine covariate"),
            HR, lower, upper)
figS_primary_sens <- make_forest(fig3a_df, "Primary and sensitivity analyses",
                                 c(0.68, 1.20), breaks = seq(0.7, 1.1, 0.1))

fig4_group_levels <- c(
  "Lactate at time zero",
  "Norepinephrine dose at time zero",
  "Baseline steroid",
  "Invasive mechanical ventilation",
  "SOFA score at time zero"
)

fig4_df <- bind_rows(
  subgroup_within %>% filter(stratum_var == "lact_stratum") %>%
    mutate(group = "Lactate at time zero",
           stratum = recode(.stratum, "<2 mmol/L" = "<2 mmol/L",
                            "2-<4 mmol/L" = "2-<4 mmol/L", ">=4 mmol/L" = ">=4 mmol/L")),
  subgroup_within %>% filter(stratum_var == "ne_high_035") %>%
    mutate(group = "Norepinephrine dose at time zero",
           stratum = if_else(.stratum == "1", ">=0.35 \u00b5g/kg/min", "<0.35 \u00b5g/kg/min")),
  subgroup_within %>% filter(stratum_var == "steroid_pre") %>%
    mutate(group = "Baseline steroid",
           stratum = if_else(.stratum == "1", "Yes", "No")),
  subgroup_within %>% filter(stratum_var == "mech_vent0_inv") %>%
    mutate(group = "Invasive mechanical ventilation",
           stratum = if_else(.stratum == "1", "Yes", "No")),
  subgroup_within %>% filter(stratum_var == "sofa_ge14") %>%
    mutate(group = "SOFA score at time zero",
           stratum = if_else(.stratum == "1", ">=14", "<14"))
) %>%
  mutate(
    group = factor(group, levels = fig4_group_levels),
    stratum_order = case_when(
      group == "Lactate at time zero" & stratum == "<2 mmol/L" ~ 1,
      group == "Lactate at time zero" & stratum == "2-<4 mmol/L" ~ 2,
      group == "Lactate at time zero" & stratum == ">=4 mmol/L" ~ 3,
      group == "Norepinephrine dose at time zero" & stratum == "<0.35 \u00b5g/kg/min" ~ 1,
      group == "Norepinephrine dose at time zero" & stratum == ">=0.35 \u00b5g/kg/min" ~ 2,
      group == "Baseline steroid" & stratum == "No" ~ 1,
      group == "Baseline steroid" & stratum == "Yes" ~ 2,
      group == "Invasive mechanical ventilation" & stratum == "No" ~ 1,
      group == "Invasive mechanical ventilation" & stratum == "Yes" ~ 2,
      group == "SOFA score at time zero" & stratum == "<14" ~ 1,
      group == "SOFA score at time zero" & stratum == ">=14" ~ 2,
      TRUE ~ 99
    )
  ) %>%
  arrange(group, stratum_order) %>%
  left_join(subgroup_interactions %>% select(sub_var, interaction_p),
            by = c("stratum_var" = "sub_var")) %>%
  group_by(group) %>%
  mutate(p_label = if_else(row_number() == 1 & !is.na(interaction_p),
                           paste0("P int = ", fmt_p(interaction_p)),
                           "")) %>%
  ungroup() %>%
  mutate(label = stratum,
         text = fmt_hr(HR, lower, upper)) %>%
  select(-stratum_order)

fig4_p_labels <- fig4_df %>%
  filter(p_label != "") %>%
  transmute(group = as.character(group), p_label)

fig4_rows <- map_dfr(fig4_group_levels, function(g) {
  bind_rows(
    tibble(group = g, label = g, is_header = TRUE, HR = NA_real_, lower = NA_real_,
           upper = NA_real_, text = NA_character_),
    fig4_df %>%
      filter(as.character(group) == g) %>%
      transmute(group = as.character(group), label, is_header = FALSE, HR, lower, upper, text)
  )
}) %>%
  left_join(fig4_p_labels, by = "group") %>%
  mutate(p_label = if_else(is_header, replace_na(p_label, ""), "")) %>%
  mutate(row = rev(row_number()))

fig4_plot_df <- fig4_rows %>% filter(!is_header)
fig4_header_df <- fig4_rows %>% filter(is_header)
fig4_ylim <- c(0.5, nrow(fig4_rows) + 0.5)

fig4_labels <- ggplot(fig4_rows, aes(y = row, label = label)) +
  geom_text(data = fig4_rows %>% filter(is_header),
            aes(x = 0.00), hjust = 0, size = 4.7, family = FONT_FAMILY,
            fontface = "bold", color = COL_BLACK) +
  geom_text(data = fig4_rows %>% filter(!is_header),
            aes(x = 0.08), hjust = 0, size = 4.3, family = FONT_FAMILY,
            fontface = "bold", color = "#4A4A4A") +
  scale_y_continuous(limits = fig4_ylim, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.margin = margin(18, 0, 52, 18))

fig4_forest <- ggplot(fig4_plot_df, aes(y = row, x = HR)) +
  geom_hline(aes(yintercept = row), color = "#E6E6E6", linewidth = 0.55) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.8, color = "#555555") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                width = 0.12, linewidth = 0.95, color = COL_BLACK) +
  geom_point(size = 4.2, color = COL_BLACK) +
  scale_y_continuous(breaks = NULL, labels = NULL,
                     limits = fig4_ylim, expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.55, 1.22), breaks = seq(0.6, 1.2, 0.2)) +
  labs(x = "Hazard ratio (ultra-early vs early)", y = NULL) +
  theme_icm(18) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 22, margin = margin(t = 12)),
        plot.margin = margin(18, 0, 18, 0))

fig4_values <- ggplot(fig4_plot_df, aes(x = 0, y = row, label = text)) +
  geom_text(hjust = 0, size = 4.35, family = FONT_FAMILY, fontface = "bold",
            color = COL_BLACK) +
  scale_y_continuous(limits = fig4_ylim, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.margin = margin(18, 18, 52, 8))

fig4_pint <- ggplot(fig4_rows, aes(x = 0, y = row, label = p_label)) +
  geom_text(data = fig4_header_df,
            hjust = 0, size = 3.75, family = FONT_FAMILY, fontface = "bold",
            color = COL_BLACK) +
  scale_y_continuous(limits = fig4_ylim, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.margin = margin(18, 18, 52, 4))

fig4 <- fig4_labels + fig4_forest + fig4_values + fig4_pint +
  plot_layout(widths = c(0.32, 0.43, 0.16, 0.09))

vp_timing_bins <- pt_ne %>%
  filter(!is.na(tvp_rel), tvp_rel <= EARLY_END) %>%
  mutate(
    timing_bin = cut(
      tvp_rel,
      breaks = c(-Inf, 1, 2, 3, 4, 5, 6),
      labels = c("0-1", ">1-2", ">2-3", ">3-4", ">4-5", ">5-6"),
      right = TRUE
    )
  ) %>%
  count(timing_bin, name = "n") %>%
  mutate(pct = n / sum(n))

figS1 <- ggplot(vp_timing_bins, aes(x = timing_bin, y = n)) +
  geom_col(fill = "white", color = COL_BLACK, linewidth = 0.8, width = 0.75) +
  geom_text(aes(label = paste0(fmt_n(n), "\n", sprintf("%.1f%%", 100 * pct))),
            vjust = -0.15, size = 5.0, family = FONT_FAMILY, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Hours from NE time zero to vasopressin initiation",
       y = "Observed patients, n") +
  theme_icm(18) +
  theme(axis.title.x = element_text(size = 22, margin = margin(t = 12)),
        axis.title.y = element_text(size = 22, margin = margin(r = 12)))

figS2 <- ccw_tbl_main %>%
  mutate(arm = factor(arm, levels = c("ultra", "early"),
                      labels = c("0-3 h", ">3-6 h")),
         w_cap = pmax(P_CLIP, pmin(weight, caps_main$cap99)))
figS2_label <- figS2 %>%
  group_by(arm) %>%
  summarise(
    med = median(w_cap, na.rm = TRUE),
    q1 = quantile(w_cap, 0.25, na.rm = TRUE),
    q3 = quantile(w_cap, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("Median %.2f (IQR %.2f-%.2f)", med, q1, q3))
figS2 <- ggplot(figS2, aes(y = arm, x = w_cap)) +
  geom_boxplot(fill = "white", color = COL_BLACK, linewidth = 0.9,
               outlier.shape = 16, outlier.size = 1.8, outlier.alpha = 0.40) +
  geom_text(data = figS2_label, aes(x = 6.2, y = arm, label = label),
            inherit.aes = FALSE, hjust = 0, size = 5.0,
            family = FONT_FAMILY, fontface = "bold") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20),
                labels = c("0.5", "1", "2", "5", "10", "20")) +
  expand_limits(x = 25) +
  labs(x = "Stabilized IPC weight after 99th-percentile truncation",
       y = NULL) +
  theme_icm(18) +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.text.y = element_text(size = 20, face = "bold"))

uncensored_counts <- w_grace_main %>%
  count(arm, t, name = "clone_intervals") %>%
  mutate(arm = factor(arm, levels = c("ultra", "early"),
                      labels = c("0-3 h", ">3-6 h")))

figS3 <- ggplot(uncensored_counts, aes(x = t, y = clone_intervals, linetype = arm)) +
  geom_line(linewidth = 1.25, color = COL_BLACK) +
  geom_point(size = 3.4, color = COL_BLACK) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_x_continuous(breaks = 0:(EARLY_END - 1)) +
  scale_y_continuous(labels = comma) +
  labs(x = "Hours after NE time zero",
       y = "Uncensored clone intervals, n", linetype = NULL) +
  theme_icm(18) +
  theme(legend.position = c(0.82, 0.82),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)))

message("Writing editable figures...")
ggsave(file.path(FIG_DIR, "figure1_flow_ne_main.pdf"), fig1, width = 11, height = 7.8, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG_DIR, "figure2_cuminc_ne_main.pdf"), fig2, width = 11, height = 7.8, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG_DIR, "figure3_outcomes_ne_main.pdf"), fig3, width = 11, height = 7.8, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG_DIR, "figure3_primary_sensitivity_ne_main.pdf"), fig3, width = 11, height = 7.8, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG_DIR, "figure4_subgroup_ne_main.pdf"), fig4, width = 11, height = 8.5, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG_DIR, "supp_figure_primary_sensitivity_ne_main.pdf"), figS_primary_sens, width = 11, height = 7.5, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG_DIR, "supp_figure1_vp_timing_ne_main.pdf"), figS1, width = 11, height = 7.5, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG_DIR, "supp_figure2_ipcw_weights_ne_main.pdf"), figS2, width = 11, height = 7.5, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG_DIR, "supp_figure3_uncensored_clones_ne_main.pdf"), figS3, width = 11, height = 7.5, device = cairo_pdf, bg = "white")
ggsave(file.path(QA_DIR, "figure1_flow_ne_main.png"), fig1, width = 11, height = 7.8, dpi = 220, bg = "white")
ggsave(file.path(QA_DIR, "figure2_cuminc_ne_main.png"), fig2, width = 11, height = 7.8, dpi = 220, bg = "white")
ggsave(file.path(QA_DIR, "figure3_outcomes_ne_main.png"), fig3, width = 11, height = 7.8, dpi = 220, bg = "white")
ggsave(file.path(QA_DIR, "figure3_primary_sensitivity_ne_main.png"), fig3, width = 11, height = 7.8, dpi = 220, bg = "white")
ggsave(file.path(QA_DIR, "figure4_subgroup_ne_main.png"), fig4, width = 11, height = 8.5, dpi = 220, bg = "white")
ggsave(file.path(QA_DIR, "supp_figure_primary_sensitivity_ne_main.png"), figS_primary_sens, width = 11, height = 7.5, dpi = 220, bg = "white")
ggsave(file.path(QA_DIR, "supp_figure1_vp_timing_ne_main.png"), figS1, width = 11, height = 7.5, dpi = 220, bg = "white")
ggsave(file.path(QA_DIR, "supp_figure2_ipcw_weights_ne_main.png"), figS2, width = 11, height = 7.5, dpi = 220, bg = "white")
ggsave(file.path(QA_DIR, "supp_figure3_uncensored_clones_ne_main.png"), figS3, width = 11, height = 7.5, dpi = 220, bg = "white")

ppt <- read_pptx()
add_figure_slide <- function(ppt, title, gg, w = 9.5, h = 6.9) {
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
  ph_with(ppt, title, location = ph_location(left = 0.25, top = 0.12, width = 9.5, height = 0.35),
          gp = fp_text(font.family = FONT_FAMILY, font.size = 12, bold = TRUE))
  ph_with(ppt, dml(ggobj = gg), location = ph_location(left = 0.25, top = 0.45, width = w, height = h))
}
ppt <- add_figure_slide(ppt, "Figure 1. Study flow diagram", fig1)
ppt <- add_figure_slide(ppt, "Figure 2. Weighted cumulative incidence curves for 28-day mortality", fig2)
ppt <- add_figure_slide(ppt, "Figure 3. Primary outcome sensitivity analyses", fig3, h = 6.9)
ppt <- add_figure_slide(ppt, "Figure 4. Subgroup analyses for 28-day mortality", fig4, h = 6.9)
fig_pptx <- file.path(OUT_DIR, "VP_NE_main_figures_editable.pptx")
print(ppt, target = fig_pptx)

supp_ppt <- read_pptx()
supp_ppt <- add_figure_slide(supp_ppt, "Supplementary Figure. Primary and sensitivity analyses", figS_primary_sens)
supp_ppt <- add_figure_slide(supp_ppt, "Supplementary Figure 1. Observed timing of vasopressin initiation", figS1)
supp_ppt <- add_figure_slide(supp_ppt, "Supplementary Figure 2. IPCW distribution", figS2)
supp_ppt <- add_figure_slide(supp_ppt, "Supplementary Figure 3. Uncensored clone intervals", figS3)
supp_fig_pptx <- file.path(OUT_DIR, "VP_NE_main_supplementary_figures_editable.pptx")
print(supp_ppt, target = supp_fig_pptx)

smd_cont <- function(x, g) {
  x1 <- x[g == levels(g)[1]]
  x2 <- x[g == levels(g)[2]]
  s <- sqrt((var(x1, na.rm = TRUE) + var(x2, na.rm = TRUE)) / 2)
  if (is.na(s) || s == 0) return(NA_real_)
  abs((mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)) / s)
}
smd_bin <- function(x, g) {
  p1 <- mean(x[g == levels(g)[1]], na.rm = TRUE)
  p2 <- mean(x[g == levels(g)[2]], na.rm = TRUE)
  s <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  if (is.na(s) || s == 0) return(NA_real_)
  abs((p1 - p2) / s)
}
smd_cont_positive <- function(x, g) {
  x1 <- x[g == levels(g)[1]]
  x2 <- x[g == levels(g)[2]]
  x1 <- x1[!is.na(x1) & x1 > 0]
  x2 <- x2[!is.na(x2) & x2 > 0]
  s <- sqrt((var(x1, na.rm = TRUE) + var(x2, na.rm = TRUE)) / 2)
  if (is.na(s) || s == 0) return(NA_real_)
  abs((mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)) / s)
}

table1_specs <- tribble(
  ~section, ~var, ~label, ~type, ~digits,
  "Demographics", "age", "Age, yr", "cont", 1,
  "Demographics", "male", "Male sex", "bin", NA_real_,
  "Demographics", "bmi", "BMI, kg/m2", "cont", 1,
  "Demographics", "wt0", "Body weight, kg", "cont", 1,
  "Hemodynamics and severity at time zero", "map0", "MAP, mmHg", "cont", 1,
  "Hemodynamics and severity at time zero", "lact0", "Lactate, mmol/L", "cont", 1,
  "Hemodynamics and severity at time zero", "sofa0", "SOFA score", "cont", 1,
  "Hemodynamics and severity at time zero", "crea0", "Creatinine, mg/dL", "cont", 2,
  "Hemodynamics and severity at time zero", "ne_rate_0", "Norepinephrine dose, mcg/kg/min", "cont", 2,
  "Hemodynamics and severity at time zero", "nee_0", "NE-equivalent dose, mcg/kg/min", "cont", 2,
  "Hemodynamics and severity at time zero", "mech_vent0_inv", "Invasive mechanical ventilation", "bin", NA_real_,
  "Co-interventions before time zero", "steroid_pre", "Systemic corticosteroid use before time zero", "bin", NA_real_,
  "Past medical history", "cc_chronic_cardiac", "Chronic cardiac disease", "bin", NA_real_,
  "Past medical history", "cc_hypertension", "Hypertension", "bin", NA_real_,
  "Past medical history", "cc_ckd", "Chronic kidney disease", "bin", NA_real_,
  "Past medical history", "cc_diabetes", "Diabetes", "bin", NA_real_,
  "Past medical history", "cc_malignancy", "Malignancy", "bin", NA_real_,
  "Past medical history", "cc_hiv", "HIV/AIDS", "bin", NA_real_,
  "Past medical history", "cc_dementia", "Dementia", "bin", NA_real_,
  "Fluid before time zero", "fluid_in_m24_t0", "Fluid input in 24 h before time zero, mL", "cont0", 0,
  "Fluid before time zero", "fluid_in_m6_t0", "Fluid input in 6 h before time zero, mL", "cont0", 0,
  "Fluid before time zero", "fluid_bal_m24_t0", "Fluid balance in 24 h before time zero, mL", "cont0", 0,
  "Fluid before time zero", "fluid_bal_m6_t0", "Fluid balance in 6 h before time zero, mL", "cont0", 0
)

pt_t1 <- pt_ne %>%
  mutate(male = sex %in% c("Male", "M", "male", "1"),
         obs2 = factor(obs_vp_group, levels = c("0-3 h", ">3-6 h"))) %>%
  filter(!is.na(obs2) | TRUE)

n_table1_overall <- nrow(pt_t1)
n_table1_ultra <- sum(pt_t1$obs_vp_group == "0-3 h", na.rm = TRUE)
n_table1_early <- sum(pt_t1$obs_vp_group == ">3-6 h", na.rm = TRUE)
n_table1_not_early <- sum(pt_t1$obs_vp_group == ">6 h", na.rm = TRUE)
n_table1_non_vp <- sum(pt_t1$obs_vp_group == "Never VP", na.rm = TRUE)
table1_col_overall <- paste0("Overall\n(N = ", fmt_n(n_table1_overall), ")")
table1_col_ultra <- paste0("Ultra-early\n0-3 h\n(n = ", fmt_n(n_table1_ultra), ")")
table1_col_early <- paste0("Early\n>3-6 h\n(n = ", fmt_n(n_table1_early), ")")
table1_col_not_early <- paste0("Not early\n>6 h\n(n = ", fmt_n(n_table1_not_early), ")")
table1_col_non_vp <- paste0("Non-VP use\n(n = ", fmt_n(n_table1_non_vp), ")")

summ_one <- function(df, var, type, digits) {
  x <- df[[var]]
  if (type == "bin") return(fmt_count_pct(as.logical(x)))
  if (type == "cont0") {
    x <- x[!is.na(x) & x > 0]
    if (!length(x)) return("")
    return(fmt_med_iqr(x, digits))
  }
  fmt_med_iqr(as.numeric(x), digits)
}

table1_body <- map_dfr(unique(table1_specs$section), function(sec) {
  specs <- table1_specs %>% filter(section == sec)
  bind_rows(
    tibble(Characteristic = sec, Overall = "", `Ultra-early 0-3 h` = "",
           `Early >3-6 h` = "", SMD = "", .section = TRUE),
    pmap_dfr(specs, function(section, var, label, type, digits) {
      gdata <- pt_t1 %>% filter(obs2 %in% c("0-3 h", ">3-6 h"))
      g <- droplevels(gdata$obs2)
      smd <- if (length(levels(g)) < 2) {
        NA_real_
      } else if (type == "bin") {
        smd_bin(as.numeric(as.logical(gdata[[var]])), g)
      } else if (type == "cont0") {
        smd_cont_positive(as.numeric(gdata[[var]]), g)
      } else {
        smd_cont(as.numeric(gdata[[var]]), g)
      }
      tibble(
        Characteristic = label,
        Overall = summ_one(pt_t1, var, type, digits),
        `Ultra-early 0-3 h` = summ_one(pt_t1 %>% filter(obs_vp_group == "0-3 h"), var, type, digits),
        `Early >3-6 h` = summ_one(pt_t1 %>% filter(obs_vp_group == ">3-6 h"), var, type, digits),
        SMD = ifelse(is.na(smd), "", sprintf("%.2f", smd)),
        .section = FALSE
      )
    })
  )
}) %>%
  bind_rows(
    tibble(
      Characteristic = "Patients, n",
      Overall = fmt_n(n_table1_overall),
      `Ultra-early 0-3 h` = fmt_n(n_table1_ultra),
      `Early >3-6 h` = fmt_n(n_table1_early),
      SMD = "",
      .section = FALSE
    ),
    .
  ) %>%
  rename(
    !!table1_col_overall := Overall,
    !!table1_col_ultra := `Ultra-early 0-3 h`,
    !!table1_col_early := `Early >3-6 h`
  )

table1_supp_body <- map_dfr(unique(table1_specs$section), function(sec) {
  specs <- table1_specs %>% filter(section == sec)
  bind_rows(
    tibble(Characteristic = sec, Overall = "", `Ultra-early 0-3 h` = "",
           `Early >3-6 h` = "", `Not early >6 h` = "",
           `Non-VP use` = "", .section = TRUE),
    pmap_dfr(specs, function(section, var, label, type, digits) {
      tibble(
        Characteristic = label,
        Overall = summ_one(pt_t1, var, type, digits),
        `Ultra-early 0-3 h` = summ_one(pt_t1 %>% filter(obs_vp_group == "0-3 h"), var, type, digits),
        `Early >3-6 h` = summ_one(pt_t1 %>% filter(obs_vp_group == ">3-6 h"), var, type, digits),
        `Not early >6 h` = summ_one(pt_t1 %>% filter(obs_vp_group == ">6 h"), var, type, digits),
        `Non-VP use` = summ_one(pt_t1 %>% filter(obs_vp_group == "Never VP"), var, type, digits),
        .section = FALSE
      )
    })
  )
}) %>%
  bind_rows(
    tibble(
      Characteristic = "Patients, n",
      Overall = fmt_n(n_table1_overall),
      `Ultra-early 0-3 h` = fmt_n(n_table1_ultra),
      `Early >3-6 h` = fmt_n(n_table1_early),
      `Not early >6 h` = fmt_n(n_table1_not_early),
      `Non-VP use` = fmt_n(n_table1_non_vp),
      .section = FALSE
    ),
    .
  ) %>%
  rename(
    !!table1_col_overall := Overall,
    !!table1_col_ultra := `Ultra-early 0-3 h`,
    !!table1_col_early := `Early >3-6 h`,
    !!table1_col_not_early := `Not early >6 h`,
    !!table1_col_non_vp := `Non-VP use`
  )

make_ft <- function(dat, title = NULL, font_size = 9) {
  ft <- dat %>% select(-any_of(".section")) %>% flextable()
  ft <- theme_vanilla(ft)
  ft <- font(ft, fontname = FONT_FAMILY, part = "all")
  ft <- fontsize(ft, size = font_size, part = "all")
  ft <- fontsize(ft, size = font_size + 1, part = "header")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "header")
  ft <- align(ft, j = 1, align = "left", part = "body")
  if (".section" %in% names(dat)) {
    sec_rows <- which(dat$.section)
    ft <- bold(ft, i = sec_rows, part = "body")
    ft <- bg(ft, i = sec_rows, bg = "#EFEFEF", part = "body")
  }
  ft <- padding(ft, padding.top = 3, padding.bottom = 3, padding.left = 4, padding.right = 4)
  ft <- autofit(ft)
  ft <- set_table_properties(ft, layout = "autofit", width = 1)
  if (!is.null(title)) {
    ft <- add_header_lines(ft, values = title)
    ft <- bold(ft, i = 1, part = "header")
    ft <- fontsize(ft, i = 1, size = font_size + 2, part = "header")
  }
  ft
}

make_ft_table1 <- function(dat, title = NULL, font_size = 8.2) {
  body_dat <- dat %>% select(-any_of(".section"))
  ft <- flextable(body_dat)
  black_thin <- fp_border(color = "black", width = 0.5)
  black_med <- fp_border(color = "black", width = 0.9)
  ft <- border_remove(ft)
  ft <- font(ft, fontname = FONT_FAMILY, part = "all")
  ft <- fontsize(ft, size = font_size, part = "all")
  ft <- fontsize(ft, size = font_size + 0.4, part = "header")
  ft <- color(ft, color = COL_BLACK, part = "all")
  ft <- bg(ft, bg = "white", part = "all")
  ft <- bold(ft, part = "header")
  ft <- bold(ft, j = 1, part = "body")
  ft <- align(ft, j = 1, align = "left", part = "all")
  ft <- align(ft, j = 2:ncol(body_dat), align = "center", part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  ft <- padding(ft, padding.top = 3, padding.bottom = 3,
                padding.left = 4, padding.right = 4, part = "all")
  ft <- hline_top(ft, border = black_med, part = "header")
  ft <- hline_bottom(ft, border = black_med, part = "header")
  ft <- hline_bottom(ft, border = black_med, part = "body")
  if (".section" %in% names(dat)) {
    sec_rows <- which(dat$.section)
    ft <- bold(ft, i = sec_rows, part = "body")
    ft <- hline(ft, i = sec_rows, border = black_thin, part = "body")
  }
  ft <- width(ft, j = 1, width = 2.35)
  ft <- width(ft, j = 2:4, width = 1.15)
  ft <- width(ft, j = 5, width = 0.50)
  ft <- set_table_properties(ft, layout = "fixed", width = 1)
  ft <- add_footer_lines(
    ft,
    values = "Overall includes patients with ultra-early VP initiation (0-3 h), early VP initiation (>3-6 h), not-early VP initiation (>6 h), and non-VP use. Values are median (IQR) or n (%). SMD compares the observed ultra-early and early VP initiation groups."
  )
  ft <- fontsize(ft, size = font_size - 0.4, part = "footer")
  ft <- italic(ft, part = "footer")
  if (!is.null(title)) {
    ft <- add_header_lines(ft, values = title)
    ft <- bold(ft, i = 1, part = "header")
    ft <- fontsize(ft, i = 1, size = font_size + 1.4, part = "header")
    ft <- align(ft, i = 1, align = "left", part = "header")
  }
  ft
}

make_ft_table1_supp <- function(dat, title = NULL, font_size = 7.4) {
  body_dat <- dat %>% select(-any_of(".section"))
  ft <- flextable(body_dat)
  black_thin <- fp_border(color = "black", width = 0.5)
  black_med <- fp_border(color = "black", width = 0.9)
  ft <- border_remove(ft)
  ft <- font(ft, fontname = FONT_FAMILY, part = "all")
  ft <- fontsize(ft, size = font_size, part = "all")
  ft <- fontsize(ft, size = font_size + 0.3, part = "header")
  ft <- color(ft, color = COL_BLACK, part = "all")
  ft <- bg(ft, bg = "white", part = "all")
  ft <- bold(ft, part = "header")
  ft <- bold(ft, j = 1, part = "body")
  ft <- align(ft, j = 1, align = "left", part = "all")
  ft <- align(ft, j = 2:ncol(body_dat), align = "center", part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  ft <- padding(ft, padding.top = 2, padding.bottom = 2,
                padding.left = 3, padding.right = 3, part = "all")
  ft <- hline_top(ft, border = black_med, part = "header")
  ft <- hline_bottom(ft, border = black_med, part = "header")
  ft <- hline_bottom(ft, border = black_med, part = "body")
  if (".section" %in% names(dat)) {
    sec_rows <- which(dat$.section)
    ft <- bold(ft, i = sec_rows, part = "body")
    ft <- hline(ft, i = sec_rows, border = black_thin, part = "body")
  }
  ft <- width(ft, j = 1, width = 2.15)
  ft <- width(ft, j = 2:ncol(body_dat), width = 0.95)
  ft <- set_table_properties(ft, layout = "fixed", width = 1)
  ft <- add_footer_lines(
    ft,
    values = "Values are median (IQR) or n (%). Not early indicates VP initiation after 6 h; non-VP use indicates no observed VP initiation. SMD was not calculated for this supplementary descriptive table."
  )
  ft <- fontsize(ft, size = font_size - 0.3, part = "footer")
  ft <- italic(ft, part = "footer")
  if (!is.null(title)) {
    ft <- add_header_lines(ft, values = title)
    ft <- bold(ft, i = 1, part = "header")
    ft <- fontsize(ft, i = 1, size = font_size + 1.2, part = "header")
    ft <- align(ft, i = 1, align = "left", part = "header")
  }
  ft
}

ft_table1 <- make_ft(
  table1_body,
  paste0("Table 1. Baseline characteristics at norepinephrine time zero (N = ",
         fmt_n(nrow(pt_ne)), ")"),
  font_size = 8.5
)

ft_table1 <- make_ft_table1(
  table1_body,
  "Table 1. Baseline characteristics at norepinephrine time zero",
  font_size = 8.1
)

ft_table1_supp <- make_ft_table1_supp(
  table1_supp_body,
  "Supplementary Table. Baseline characteristics by observed vasopressin timing",
  font_size = 7.3
)

sens_table <- primary_sens %>%
  transmute(Analysis = label, `HR (95% CI)` = HR_95CI,
            `Clones, n` = fmt_n(n_clone), `Events, n` = fmt_n(n_event))
ft_sens <- make_ft(sens_table, "Primary and key sensitivity analyses", 9)

secondary_table <- secondary_finegray %>%
  transmute(Outcome = label, `HR (95% CI)` = HR_95CI,
            `Clones, n` = fmt_n(n_clone), `Events, n` = fmt_n(n_event))
ft_secondary <- make_ft(secondary_table, "Secondary outcomes", 9)

first_vaso_table <- pt_ne %>%
  mutate(src = recode(src, miiv = "MIMIC-IV", eicu = "eICU-CRD"),
         first_vaso = replace_na(as.character(first_vaso), "missing/unknown")) %>%
  count(src, first_vaso, name = "n") %>%
  group_by(src) %>%
  mutate(`%` = 100 * n / sum(n)) %>%
  ungroup() %>%
  arrange(src, desc(n)) %>%
  transmute(Database = src, `First vasopressor` = first_vaso, n = fmt_n(n), `%` = sprintf("%.1f", `%`))
ft_first_vaso <- make_ft(first_vaso_table, "First vasopressor distribution by database", 9)

vp_dose_table <- pt_ne %>%
  filter(!is.na(tvp_rel), tvp_rel <= 6) %>%
  mutate(timing = if_else(tvp_rel <= 3, "0-3 h", ">3-6 h")) %>%
  group_by(timing) %>%
  summarise(
    n = n(),
    `Initial dose, U/min` = fmt_med_iqr(vp_dose_init, 3),
    `Maximum dose within 24 h, U/min` = fmt_med_iqr(vp_dose_max_24h, 3),
    `Cumulative dose within 24 h, U` = fmt_med_iqr(vp_cum_24h, 2),
    `Initial dose <=0.03 U/min` = fmt_count_pct(coalesce(vp_init_eq003, FALSE) | coalesce(vp_init_lt003, FALSE)),
    .groups = "drop"
  ) %>%
  mutate(n = fmt_n(n))
ft_vp_dose <- make_ft(vp_dose_table, "Vasopressin dosing by observed timing", 9)

vp_interrupt_table <- pt_ne %>%
  filter(!is.na(tvp_rel)) %>%
  mutate(timing = case_when(tvp_rel <= 3 ~ "0-3 h",
                            tvp_rel > 3 & tvp_rel <= 6 ~ ">3-6 h",
                            tvp_rel > 6 ~ ">6 h")) %>%
  group_by(timing) %>%
  summarise(
    n = n(),
    `3-h interruption` = fmt_count_pct(coalesce(interrupt_3h, FALSE)),
    `Sustained 6-h interruption` = fmt_count_pct(coalesce(sustained_interrupt_6h, FALSE)),
    `3-h permanent discontinuation` = fmt_count_pct(coalesce(permanent_disc_3h, FALSE)),
    .groups = "drop"
  ) %>%
  mutate(n = fmt_n(n))
ft_vp_interrupt <- make_ft(vp_interrupt_table, "Early vasopressin interruption", 9)

new_inotrope_table <- pt_ne %>%
  filter(!is.na(tvp_rel), tvp_rel <= 6, !coalesce(dobu_at_vp, FALSE)) %>%
  mutate(timing = if_else(tvp_rel <= 3, "0-3 h", ">3-6 h")) %>%
  group_by(timing) %>%
  summarise(
    n = n(),
    `New dobutamine within 3 h` = fmt_count_pct(coalesce(new_dobu_3h, FALSE)),
    `New dobutamine within 6 h` = fmt_count_pct(coalesce(new_dobu_6h, FALSE)),
    `New dobutamine within 24 h` = fmt_count_pct(coalesce(new_dobu_24h, FALSE)),
    `New epinephrine within 24 h` = fmt_count_pct(coalesce(new_epi_24h, FALSE)),
    .groups = "drop"
  ) %>%
  mutate(n = fmt_n(n))
ft_new_inotrope <- make_ft(new_inotrope_table, "New inotrope or epinephrine initiation after VP", 9)

no_vp_table <- pt_ne %>%
  group_by(obs_vp_group) %>%
  summarise(
    n = n(),
    `Age, yr` = fmt_med_iqr(age, 1),
    `SOFA score` = fmt_med_iqr(sofa0, 1),
    `Lactate, mmol/L` = fmt_med_iqr(lact0, 1),
    `Norepinephrine dose, mcg/kg/min` = fmt_med_iqr(ne_rate_0, 2),
    `Invasive ventilation` = fmt_count_pct(as.logical(mech_vent0_inv)),
    `Crude 28-d mortality` = fmt_count_pct(death28),
    .groups = "drop"
  ) %>%
  rename(`Observed VP timing` = obs_vp_group) %>%
  mutate(n = fmt_n(n))
ft_no_vp <- make_ft(no_vp_table, "Observed VP timing groups, including never VP", 9)

safety_table <- pt_ne %>%
  group_by(obs_vp_group) %>%
  summarise(
    n = n(),
    DVT = fmt_count_pct(coalesce(sf_dvt, FALSE)),
    PE = fmt_count_pct(coalesce(sf_pe, FALSE)),
    `Mesenteric ischemia` = fmt_count_pct(coalesce(sf_mesenteric, FALSE)),
    `Digital/peripheral ischemia` = fmt_count_pct(coalesce(sf_digital, FALSE)),
    .groups = "drop"
  ) %>%
  rename(`Observed VP timing` = obs_vp_group) %>%
  mutate(n = fmt_n(n))
ft_safety <- make_ft(safety_table, "Exploratory diagnosis-based safety outcomes", 9)

fluid_vars <- c(
  fluid_in_m24_t0 = "Fluid input in 24 h before time zero",
  fluid_in_m6_t0 = "Fluid input in 6 h before time zero",
  fluid_bal_m24_t0 = "Fluid balance in 24 h before time zero",
  fluid_bal_m6_t0 = "Fluid balance in 6 h before time zero",
  fluid_in_m24_vp = "Fluid input in 24 h before VP",
  fluid_in_m6_vp = "Fluid input in 6 h before VP",
  fluid_in_p6_vp = "Fluid input in 6 h after VP",
  fluid_bal_m24_vp = "Fluid balance in 24 h before VP",
  fluid_bal_m6_vp = "Fluid balance in 6 h before VP",
  fluid_bal_p6_vp = "Fluid balance in 6 h after VP"
)
fluid_summary <- function(df, var) {
  x <- df[[var]]
  x <- x[!is.na(x) & x > 0]
  if (!length(x)) return("")
  sprintf("%s [n=%s]", fmt_med_iqr(x, 0), fmt_n(length(x)))
}
fluid_table <- crossing(
  timing = c("0-3 h", ">3-6 h"),
  variable = names(fluid_vars)
) %>%
  rowwise() %>%
  mutate(
    Summary = fluid_summary(pt_ne %>% filter(obs_vp_group == timing), variable),
    Variable = unname(fluid_vars[variable])
  ) %>%
  ungroup() %>%
  select(`Observed VP timing` = timing, Variable, Summary)
ft_fluid <- make_ft(fluid_table, "Fluid input and balance around time zero and VP", 8.5)

first_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  x[1]
}
fmt_count_total <- function(n, total) {
  if (is.na(total) || total <= 0) return("")
  sprintf("%s / %s (%.1f%%)", fmt_n(n), fmt_n(total), 100 * n / total)
}

hemo_delay_base <- pt_ne %>%
  filter(!is.na(tvp_rel), tvp_rel > 3, tvp_rel <= 6) %>%
  select(uid, t0_hr, map0, ne_rate_0)

hemo_delay_panel <- ts_imputed %>%
  inner_join(hemo_delay_base %>% select(uid, t0_hr), by = "uid") %>%
  mutate(t_rel = as.integer(round(time_hr - t0_hr))) %>%
  filter(t_rel >= 0, t_rel <= 3) %>%
  group_by(uid, t_rel) %>%
  summarise(
    map = first(map),
    ne_rate = first(norepi_rate),
    epi_rate = first(epi_rate),
    dopa_rate = first(dopa_rate),
    dobu_rate = first(dobu_rate),
    phn_rate = first(phn_rate),
    mech_vent = first(mech_vent),
    .groups = "drop"
  )

hemo_delay_patient <- hemo_delay_panel %>%
  group_by(uid) %>%
  summarise(
    min_map_0_3 = if (all(is.na(map))) NA_real_ else min(map, na.rm = TRUE),
    mean_map_0_3 = if (all(is.na(map))) NA_real_ else mean(map, na.rm = TRUE),
    ne_rate_end_3 = first_nonmissing(ne_rate[t_rel == 3]),
    ne_rate_max_0_3 = if (all(is.na(ne_rate))) NA_real_ else max(ne_rate, na.rm = TRUE),
    any_epi_0_3 = any(replace_na(epi_rate, 0) > 0),
    any_dopa_0_3 = any(replace_na(dopa_rate, 0) > 0),
    any_dobu_0_3 = any(replace_na(dobu_rate, 0) > 0),
    any_phn_0_3 = any(replace_na(phn_rate, 0) > 0),
    any_inv_vent_0_3 = any(replace_na(mech_vent == "invasive", FALSE)),
    .groups = "drop"
  ) %>%
  right_join(hemo_delay_base, by = "uid") %>%
  mutate(ne_rate_delta_max = ne_rate_max_0_3 - ne_rate_0)

hemo_delay_n <- nrow(hemo_delay_base)
hemo_delay_table <- tibble(
  Characteristic = c(
    "Patients, n",
    "MAP at time zero, mmHg",
    "Minimum MAP during first 3 h, mmHg",
    "Mean MAP during first 3 h, mmHg",
    "Norepinephrine dose at time zero, mcg/kg/min",
    "Norepinephrine dose at end of 0-3 h window, mcg/kg/min",
    "Maximum norepinephrine dose during first 3 h, mcg/kg/min",
    "Increase in norepinephrine dose to maximum (0-3 h), mcg/kg/min",
    "Any epinephrine use during first 3 h",
    "Any dopamine use during first 3 h",
    "Any dobutamine use during first 3 h",
    "Any phenylephrine use during first 3 h",
    "Any invasive ventilation during first 3 h"
  ),
  Value = c(
    fmt_n(hemo_delay_n),
    fmt_med_iqr(hemo_delay_patient$map0, 1),
    fmt_med_iqr(hemo_delay_patient$min_map_0_3, 1),
    fmt_med_iqr(hemo_delay_patient$mean_map_0_3, 1),
    fmt_med_iqr(hemo_delay_patient$ne_rate_0, 2),
    fmt_med_iqr(hemo_delay_patient$ne_rate_end_3, 2),
    fmt_med_iqr(hemo_delay_patient$ne_rate_max_0_3, 2),
    fmt_med_iqr(hemo_delay_patient$ne_rate_delta_max, 2),
    fmt_count_total(sum(hemo_delay_patient$any_epi_0_3, na.rm = TRUE), hemo_delay_n),
    fmt_count_total(sum(hemo_delay_patient$any_dopa_0_3, na.rm = TRUE), hemo_delay_n),
    fmt_count_total(sum(hemo_delay_patient$any_dobu_0_3, na.rm = TRUE), hemo_delay_n),
    fmt_count_total(sum(hemo_delay_patient$any_phn_0_3, na.rm = TRUE), hemo_delay_n),
    fmt_count_total(sum(hemo_delay_patient$any_inv_vent_0_3, na.rm = TRUE), hemo_delay_n)
  )
)
ft_hemo_delay <- make_ft(hemo_delay_table, "Hemodynamic co-interventions in the observed >3-6 h group during the first 0-3 h", 9)

steroid_timing_table <- pt_ne %>%
  filter(!is.na(tvp_rel), tvp_rel <= 6) %>%
  mutate(
    timing = if_else(tvp_rel <= 3, "0-3 h", ">3-6 h"),
    steroid_cat = case_when(
      coalesce(steroid_pre_vp, FALSE) ~ "Before VP",
      !coalesce(steroid_pre_vp, FALSE) & coalesce(steroid_post_vp, FALSE) ~ "After VP",
      TRUE ~ "None"
    )
  ) %>%
  count(timing, steroid_cat, name = "n") %>%
  group_by(timing) %>%
  mutate(`%` = 100 * n / sum(n)) %>%
  ungroup() %>%
  transmute(`Observed VP timing` = timing, `Steroid timing` = steroid_cat,
            n = fmt_n(n), `%` = sprintf("%.1f", `%`))
ft_steroid <- make_ft(steroid_timing_table, "Steroid timing relative to VP", 9)

survivor_table <- pt_ne %>%
  mutate(surv28 = factor(if_else(death28, "Non-survivor", "Survivor"),
                         levels = c("Survivor", "Non-survivor")),
         male = sex %in% c("Male", "M", "male", "1")) %>%
  group_by(surv28) %>%
  summarise(
    n = n(),
    `Age, yr` = fmt_med_iqr(age, 1),
    `Male sex` = fmt_count_pct(male),
    `MAP, mmHg` = fmt_med_iqr(map0, 1),
    `Lactate, mmol/L` = fmt_med_iqr(lact0, 1),
    `SOFA score` = fmt_med_iqr(sofa0, 1),
    `Creatinine, mg/dL` = fmt_med_iqr(crea0, 2),
    `Norepinephrine dose, mcg/kg/min` = fmt_med_iqr(ne_rate_0, 2),
    `Chronic cardiac disease` = fmt_count_pct(cc_chronic_cardiac),
    `Hypertension` = fmt_count_pct(cc_hypertension),
    `Chronic kidney disease` = fmt_count_pct(cc_ckd),
    `Malignancy` = fmt_count_pct(cc_malignancy),
    .groups = "drop"
  ) %>%
  rename(`28-d status` = surv28) %>%
  mutate(n = fmt_n(n))
ft_survivor <- make_ft(survivor_table, "Survivors versus non-survivors at 28 days", 8.5)

caption_entries <- tribble(
  ~Item, ~Caption, ~Abbreviations,
  "Table 1",
  "Baseline characteristics at norepinephrine time zero. Overall includes patients with ultra-early vasopressin initiation, early vasopressin initiation, not-early vasopressin initiation, and non-vasopressin use; SMDs compare the observed ultra-early and early vasopressin initiation groups.",
  "BMI, body mass index; IQR, interquartile range; NE, norepinephrine; SMD, standardized mean difference; VP, vasopressin.",
  "Figure 1",
  "Study flow diagram for the clone-censor-weight target trial emulation. Eligible patients were cloned into the ultra-early and early vasopressin initiation strategies and artificially censored when their observed treatment deviated from the assigned strategy.",
  "CRD, Collaborative Research Database; eICU-CRD, eICU Collaborative Research Database; MIMIC-IV, Medical Information Mart for Intensive Care IV; NE, norepinephrine; VP, vasopressin.",
  "Figure 2",
  "Weighted cumulative incidence curves for 28-day mortality according to the ultra-early and early vasopressin initiation strategies. Curves were estimated after stabilized inverse-probability-of-censoring weighting with 99th-percentile truncation.",
  "CI, confidence interval; IPCW, inverse-probability-of-censoring weighting; RD, risk difference; RMST, restricted mean survival time.",
  "Figure 3",
  "Primary outcome sensitivity analysis comparing the 99th- and 95th-percentile truncation thresholds for stabilized inverse-probability-of-censoring weights.",
  "CI, confidence interval; HR, hazard ratio; IPCW, inverse-probability-of-censoring weighting.",
  "Figure 4",
  "Subgroup analyses for 28-day mortality. Hazard ratios compare the ultra-early strategy with the early strategy within each prespecified subgroup; P values are from weighted Cox models including treatment-by-subgroup interaction terms.",
  "CI, confidence interval; HR, hazard ratio; P int, P value for interaction; SOFA, Sequential Organ Failure Assessment.",
  "Supplementary Table 1",
  "Baseline characteristics by observed vasopressin timing group. This descriptive table separates patients into ultra-early, early, not-early, and non-vasopressin-use groups; SMDs were not calculated.",
  "BMI, body mass index; IQR, interquartile range; NE, norepinephrine; VP, vasopressin.",
  "Supplementary Table 2",
  "Primary and key sensitivity analyses for 28-day mortality.",
  "CI, confidence interval; HR, hazard ratio; IPCW, inverse-probability-of-censoring weighting; NE, norepinephrine; NEE, norepinephrine-equivalent dose.",
  "Supplementary Table 3",
  "Secondary outcomes estimated using weighted Fine-Gray models, treating death as a competing event where applicable.",
  "AKI, acute kidney injury; CI, confidence interval; CRRT, continuous renal replacement therapy; HR, hazard ratio; RRT, renal replacement therapy.",
  "Supplementary Table 4",
  "First vasopressor distribution by database at the relevant time-zero definition.",
  "eICU-CRD, eICU Collaborative Research Database; MIMIC-IV, Medical Information Mart for Intensive Care IV; NE, norepinephrine.",
  "Supplementary Table 5",
  "Vasopressin dosing among patients with observed vasopressin initiation within the 0-to-6-hour assignment window.",
  "IQR, interquartile range; VP, vasopressin.",
  "Supplementary Table 6",
  "Early vasopressin interruption after observed vasopressin initiation.",
  "VP, vasopressin.",
  "Supplementary Table 7",
  "New inotrope or epinephrine initiation after observed vasopressin initiation.",
  "VP, vasopressin.",
  "Supplementary Table 8",
  "Observed vasopressin timing groups, including patients without observed vasopressin initiation.",
  "VP, vasopressin.",
  "Supplementary Table 9",
  "Exploratory diagnosis-based safety outcomes.",
  "CI, confidence interval; HR, hazard ratio.",
  "Supplementary Table 10",
  "Fluid input and fluid balance around time zero and observed vasopressin initiation.",
  "IQR, interquartile range; VP, vasopressin.",
  "Supplementary Table 11",
  "Hemodynamic co-interventions among patients with observed early vasopressin initiation during the first 0-to-3-hour period after time zero.",
  "IQR, interquartile range; MAP, mean arterial pressure; NE, norepinephrine.",
  "Supplementary Table 12",
  "Systemic corticosteroid timing relative to observed vasopressin initiation.",
  "VP, vasopressin.",
  "Supplementary Table 13",
  "Baseline characteristics of 28-day survivors and non-survivors.",
  "IQR, interquartile range; MAP, mean arterial pressure; NE, norepinephrine; SOFA, Sequential Organ Failure Assessment.",
  "Supplementary Figure",
  "Extended primary and sensitivity analyses for 28-day mortality.",
  "CI, confidence interval; HR, hazard ratio; NE, norepinephrine; NEE, norepinephrine-equivalent dose.",
  "Supplementary Figure 1",
  "Observed timing of vasopressin initiation within 6 hours after norepinephrine time zero.",
  "NE, norepinephrine; VP, vasopressin.",
  "Supplementary Figure 2",
  "Distribution of stabilized inverse-probability-of-censoring weights after 99th-percentile truncation.",
  "IPCW, inverse-probability-of-censoring weighting; IQR, interquartile range.",
  "Supplementary Figure 3",
  "Number of uncensored clone intervals during the 0-to-6-hour assignment window.",
  "NE, norepinephrine."
)

abbreviation_entries <- tribble(
  ~Abbreviation, ~Definition,
  "AKI", "acute kidney injury",
  "BMI", "body mass index",
  "CI", "confidence interval",
  "CRD", "Collaborative Research Database",
  "CRRT", "continuous renal replacement therapy",
  "eICU-CRD", "eICU Collaborative Research Database",
  "HR", "hazard ratio",
  "ICU", "intensive care unit",
  "IPCW", "inverse-probability-of-censoring weighting",
  "IQR", "interquartile range",
  "MAP", "mean arterial pressure",
  "MIMIC-IV", "Medical Information Mart for Intensive Care IV",
  "NE", "norepinephrine",
  "NEE", "norepinephrine-equivalent dose",
  "P int", "P value for interaction",
  "RD", "risk difference",
  "RMST", "restricted mean survival time",
  "RRT", "renal replacement therapy",
  "SMD", "standardized mean difference",
  "SOFA", "Sequential Organ Failure Assessment",
  "VP", "vasopressin"
)

make_caption_ft <- function(dat, font_size = 8.2) {
  ft <- flextable(dat)
  ft <- theme_vanilla(ft)
  ft <- font(ft, fontname = FONT_FAMILY, part = "all")
  ft <- fontsize(ft, size = font_size, part = "all")
  ft <- fontsize(ft, size = font_size + 0.6, part = "header")
  ft <- bold(ft, part = "header")
  ft <- align(ft, j = 1, align = "left", part = "all")
  ft <- align(ft, j = 2:ncol(dat), align = "left", part = "all")
  ft <- valign(ft, valign = "top", part = "all")
  ft <- padding(ft, padding.top = 3, padding.bottom = 3,
                padding.left = 4, padding.right = 4, part = "all")
  ft <- width(ft, j = 1, width = 1.15)
  if (ncol(dat) == 3) {
    ft <- width(ft, j = 2, width = 4.0)
    ft <- width(ft, j = 3, width = 2.2)
  } else {
    ft <- width(ft, j = 2, width = 5.8)
  }
  ft <- set_table_properties(ft, layout = "fixed", width = 1)
  ft
}

message("Writing editable tables...")
tables_docx <- file.path(OUT_DIR, "VP_NE_main_tables_editable.docx")
save_as_docx(
  `Table 1` = ft_table1,
  `Supplementary baseline characteristics by observed VP timing` = ft_table1_supp,
  `Primary and sensitivity analyses` = ft_sens,
  `Secondary outcomes` = ft_secondary,
  `First vasopressor distribution` = ft_first_vaso,
  `Vasopressin dosing` = ft_vp_dose,
  `VP interruption` = ft_vp_interrupt,
  `New inotrope after VP` = ft_new_inotrope,
  `Observed timing including never VP` = ft_no_vp,
  `Safety outcomes` = ft_safety,
  `Fluids` = ft_fluid,
  `Hemodynamic co-interventions` = ft_hemo_delay,
  `Steroid timing` = ft_steroid,
  `Survivor versus non-survivor` = ft_survivor,
  path = tables_docx
)
save_as_docx(`Table 1` = ft_table1, path = file.path(TAB_DIR, "table1_ne_main.docx"))
table1_final_doc <- read_docx()
table1_final_doc <- body_add_flextable(table1_final_doc, ft_table1)
print(table1_final_doc, target = file.path(TAB_DIR, "table1_ne_main_final.docx"))
table1_supp_doc <- read_docx()
table1_supp_doc <- body_add_flextable(table1_supp_doc, ft_table1_supp)
print(table1_supp_doc, target = file.path(TAB_DIR, "supp_table_baseline_by_observed_vp_timing.docx"))
captions_docx <- file.path(OUT_DIR, "VP_NE_main_captions_abbreviations.docx")
captions_doc <- read_docx()
captions_doc <- body_add_par(captions_doc, "Captions and Abbreviations", style = "heading 1")
captions_doc <- body_add_par(captions_doc, "Table and figure captions", style = "heading 2")
captions_doc <- body_add_flextable(captions_doc, make_caption_ft(caption_entries, 7.6))
captions_doc <- body_add_par(captions_doc, "Abbreviations", style = "heading 2")
captions_doc <- body_add_flextable(captions_doc, make_caption_ft(abbreviation_entries, 8.5))
print(captions_doc, target = captions_docx)
save_as_docx(`Primary and sensitivity analyses` = ft_sens, path = file.path(TAB_DIR, "table_primary_sensitivity_ne_main.docx"))
save_as_docx(`Secondary outcomes` = ft_secondary, path = file.path(TAB_DIR, "table_secondary_ne_main.docx"))
save_as_docx(`Baseline characteristics by observed VP timing` = ft_table1_supp,
             `Supplementary reviewer-requested tables` = ft_first_vaso,
             `Vasopressin dosing` = ft_vp_dose,
             `VP interruption` = ft_vp_interrupt,
             `New inotrope after VP` = ft_new_inotrope,
             `Observed timing including never VP` = ft_no_vp,
             `Safety outcomes` = ft_safety,
             `Fluids` = ft_fluid,
             `Hemodynamic co-interventions` = ft_hemo_delay,
             `Steroid timing` = ft_steroid,
             `Survivor versus non-survivor` = ft_survivor,
             path = file.path(TAB_DIR, "supplementary_tables_ne_main.docx"))

message("Writing revised Online Resources draft...")
online_docx <- file.path(OUT_DIR, "VP_NE_main_online_resources_draft.docx")
or_doc <- read_docx()
or_doc <- body_add_par(or_doc, "Online Resources", style = "heading 1")
or_doc <- body_add_par(or_doc, "Ultra-early versus early adjunctive vasopressin initiation after norepinephrine escalation in septic shock", style = "Normal")

add_or_heading <- function(doc, title) {
  body_add_par(doc, title, style = "heading 2")
}
add_or_note <- function(doc, text) {
  body_add_par(doc, text, style = "Normal")
}

or_doc <- add_or_heading(or_doc, "Online Resource 1. Target trial emulation and clone-censor-weight framework")
or_doc <- add_or_note(or_doc, paste0(
  "The revised primary analysis defines time zero as the first hour at which the norepinephrine infusion rate reached ",
  "0.25 mcg/kg/min, with no prior vasopressin. Each eligible patient was cloned into the ultra-early (0-3 h) and early (>3-6 h) strategies. ",
  "Artificial censoring occurred when observed vasopressin initiation deviated from the assigned strategy. Stabilized inverse-probability-of-censoring weights were estimated using baseline and hourly time-varying covariates during the assignment window and truncated at the 99th percentile."
))

or_doc <- add_or_heading(or_doc, "Online Resource 2. Primary and key sensitivity analyses")
or_doc <- body_add_flextable(or_doc, ft_sens)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 3. Secondary outcomes")
or_doc <- body_add_flextable(or_doc, ft_secondary)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 4. First vasopressor distribution by database")
or_doc <- body_add_flextable(or_doc, ft_first_vaso)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 5. Vasopressin dosing by observed timing")
or_doc <- body_add_flextable(or_doc, ft_vp_dose)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 6. Early vasopressin interruption")
or_doc <- body_add_flextable(or_doc, ft_vp_interrupt)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 7. New inotrope or epinephrine initiation after vasopressin")
or_doc <- body_add_flextable(or_doc, ft_new_inotrope)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 8. Steroid timing relative to vasopressin")
or_doc <- body_add_flextable(or_doc, ft_steroid)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 9. Observed vasopressin timing groups, including never vasopressin")
or_doc <- add_or_note(or_doc, "This table is descriptive and is not used as a causal comparison, because never vasopressin patients are expected to differ strongly by indication and clinical trajectory.")
or_doc <- body_add_flextable(or_doc, ft_no_vp)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 10. Exploratory diagnosis-based safety outcomes")
or_doc <- body_add_flextable(or_doc, ft_safety)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 11. Fluid input and balance around time zero and vasopressin")
or_doc <- body_add_flextable(or_doc, ft_fluid)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 12. Baseline characteristics by observed vasopressin timing")
or_doc <- body_add_flextable(or_doc, ft_table1_supp)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 13. Survivors versus non-survivors at 28 days")
or_doc <- body_add_flextable(or_doc, ft_survivor)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 14. Observed timing of vasopressin initiation within 6 hours")
or_doc <- body_add_img(or_doc, src = file.path(QA_DIR, "supp_figure1_vp_timing_ne_main.png"), width = 6.4, height = 4.4)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 15. Distribution of truncated stabilized IPC weights")
or_doc <- body_add_img(or_doc, src = file.path(QA_DIR, "supp_figure2_ipcw_weights_ne_main.png"), width = 6.4, height = 4.4)
or_doc <- body_add_par(or_doc, "")

or_doc <- add_or_heading(or_doc, "Online Resource 16. Uncensored clone intervals during the 0-to-6-hour assignment window")
or_doc <- body_add_img(or_doc, src = file.path(QA_DIR, "supp_figure3_uncensored_clones_ne_main.png"), width = 6.4, height = 4.4)
print(or_doc, target = online_docx)

write_csv(table1_body %>% select(-.section), file.path(TAB_DIR, "table1_ne_main.csv"))
write_csv(table1_supp_body %>% select(-.section), file.path(TAB_DIR, "supp_table_baseline_by_observed_vp_timing.csv"))
write_csv(caption_entries, file.path(OUT_DIR, "VP_NE_main_captions.csv"))
write_csv(abbreviation_entries, file.path(OUT_DIR, "VP_NE_main_abbreviations.csv"))
write_csv(first_vaso_table, file.path(TAB_DIR, "table_first_vasopressor_ne_main.csv"))
write_csv(vp_dose_table, file.path(TAB_DIR, "table_vp_dose_ne_main.csv"))
write_csv(vp_interrupt_table, file.path(TAB_DIR, "table_vp_interruption_ne_main.csv"))
write_csv(new_inotrope_table, file.path(TAB_DIR, "table_new_inotrope_after_vp_ne_main.csv"))
write_csv(no_vp_table, file.path(TAB_DIR, "table_observed_timing_ne_main.csv"))
write_csv(safety_table, file.path(TAB_DIR, "table_safety_ne_main.csv"))
write_csv(fluid_table, file.path(TAB_DIR, "table_fluid_ne_main.csv"))
write_csv(hemo_delay_table, file.path(TAB_DIR, "table_hemodynamic_codelay_ne_main.csv"))
write_csv(steroid_timing_table, file.path(TAB_DIR, "table_steroid_timing_ne_main.csv"))
write_csv(survivor_table, file.path(TAB_DIR, "table_survivor_non_survivor_ne_main.csv"))

message("Writing requested supplementary tables and figures bundle...")
SUPP_BUNDLE_DIR <- file.path(OUT_DIR, "supplementary_tables_figures_bundle")
dir_create(SUPP_BUNDLE_DIR)

add_ft_note <- function(ft, note, font_size = 7.2) {
  ft <- add_footer_lines(ft, values = note)
  ft <- fontsize(ft, size = font_size, part = "footer")
  italic(ft, part = "footer")
}

make_timing_label <- function(x) {
  factor(
    x,
    levels = c("0-3 h", ">3-6 h", ">6 h", "Never VP"),
    labels = c("Ultra-early (0-3 h)", "Early (>3-6 h)",
               "Not early (>6 h)", "Non-VP use")
  )
}

n_table1_other <- n_table1_not_early + n_table1_non_vp
table1_col_other <- paste0("Other\n(>6 h or no VP)\n(n = ", fmt_n(n_table1_other), ")")
table1_three_body <- map_dfr(unique(table1_specs$section), function(sec) {
  specs <- table1_specs %>% filter(section == sec)
  bind_rows(
    tibble(Characteristic = sec, Overall = "", `Ultra-early 0-3 h` = "",
           `Early >3-6 h` = "", Other = "", .section = TRUE),
    pmap_dfr(specs, function(section, var, label, type, digits) {
      tibble(
        Characteristic = label,
        Overall = summ_one(pt_t1, var, type, digits),
        `Ultra-early 0-3 h` = summ_one(pt_t1 %>% filter(obs_vp_group == "0-3 h"), var, type, digits),
        `Early >3-6 h` = summ_one(pt_t1 %>% filter(obs_vp_group == ">3-6 h"), var, type, digits),
        Other = summ_one(pt_t1 %>% filter(obs_vp_group %in% c(">6 h", "Never VP")), var, type, digits),
        .section = FALSE
      )
    })
  )
}) %>%
  bind_rows(
    tibble(
      Characteristic = "Patients, n",
      Overall = fmt_n(n_table1_overall),
      `Ultra-early 0-3 h` = fmt_n(n_table1_ultra),
      `Early >3-6 h` = fmt_n(n_table1_early),
      Other = fmt_n(n_table1_other),
      .section = FALSE
    ),
    .
  ) %>%
  rename(
    !!table1_col_overall := Overall,
    !!table1_col_ultra := `Ultra-early 0-3 h`,
    !!table1_col_early := `Early >3-6 h`,
    !!table1_col_other := Other
  )
ft_baseline_three <- make_ft_table1_supp(
  table1_three_body,
  "Supplementary Table. Baseline characteristics: ultra-early, early, and other observed groups",
  font_size = 7.5
) %>%
  add_ft_note(
    "Values are median (IQR) or n (%). Other includes patients with VP initiation after 6 h and patients without observed VP initiation.",
    7.1
  )

requested_sens <- bind_rows(
  main_fit$result %>%
    mutate(Analysis = "Primary: NE time zero",
           n_patients = nrow(pt_ne),
           n_ultra = sum(!is.na(pt_ne$tvp_rel) & pt_ne$tvp_rel <= 3),
           n_early = sum(!is.na(pt_ne$tvp_rel) & pt_ne$tvp_rel > 3 & pt_ne$tvp_rel <= 6)),
  nee_fit$result %>%
    mutate(Analysis = "Sensitivity: NEE time zero",
           n_patients = nrow(pt_nee),
           n_ultra = sum(!is.na(pt_nee$tvp_rel) & pt_nee$tvp_rel <= 3),
           n_early = sum(!is.na(pt_nee$tvp_rel) & pt_nee$tvp_rel > 3 & pt_nee$tvp_rel <= 6)),
  ne_alone_fit$result %>%
    mutate(Analysis = "Sensitivity: Nad-only time zero",
           n_patients = nrow(pt_ne_alone),
           n_ultra = sum(!is.na(pt_ne_alone$tvp_rel) & pt_ne_alone$tvp_rel <= 3),
           n_early = sum(!is.na(pt_ne_alone$tvp_rel) & pt_ne_alone$tvp_rel > 3 & pt_ne_alone$tvp_rel <= 6)),
  no_dobu_t0_fit$result %>%
    mutate(Analysis = "Sensitivity: excluding dobutamine at time zero",
           n_patients = sum(!coalesce(pt_ne$dobu_at_t0, FALSE)),
           n_ultra = sum(!coalesce(pt_ne$dobu_at_t0, FALSE) & !is.na(pt_ne$tvp_rel) & pt_ne$tvp_rel <= 3),
           n_early = sum(!coalesce(pt_ne$dobu_at_t0, FALSE) & !is.na(pt_ne$tvp_rel) & pt_ne$tvp_rel > 3 & pt_ne$tvp_rel <= 6))
) %>%
  transmute(Analysis, HR, lower, upper, n_patients, n_ultra, n_early,
            `HR (95% CI)` = fmt_hr(HR, lower, upper),
            `Patients, n` = fmt_n(n_patients),
            `Observed VP 0-3 h, n` = fmt_n(n_ultra),
            `Observed VP >3-6 h, n` = fmt_n(n_early))
requested_sens_table <- requested_sens %>%
  select(Analysis, `Patients, n`, `Observed VP 0-3 h, n`,
         `Observed VP >3-6 h, n`, `HR (95% CI)`)
ft_requested_sens <- make_ft(
  requested_sens_table,
  "Supplementary Table. Requested sensitivity analyses for 28-day mortality",
  8.6
) %>%
  add_ft_note(
    "Hazard ratios compare the ultra-early strategy with the early strategy using the clone-censor-weight framework and 99th-percentile weight truncation.",
    7.2
  )

requested_sens_fig_df <- requested_sens %>%
  transmute(
    label = Analysis,
    n_label = paste0("N = ", fmt_n(n_patients),
                     "; observed 0-3 h / >3-6 h = ",
                     fmt_n(n_ultra), " / ", fmt_n(n_early)),
    HR, lower, upper
  ) %>%
  mutate(row = rev(row_number()), text = fmt_hr(HR, lower, upper))
requested_sens_ylim <- c(0.5, nrow(requested_sens_fig_df) + 0.5)
fig_requested_sens_labels <- ggplot(requested_sens_fig_df, aes(y = row)) +
  geom_text(aes(x = 0, label = label), hjust = 0, size = 4.35,
            family = FONT_FAMILY, fontface = "bold", color = COL_BLACK) +
  geom_text(aes(x = 0, y = row - 0.20, label = n_label), hjust = 0, size = 3.75,
            family = FONT_FAMILY, fontface = "bold", color = "#4A4A4A") +
  scale_y_continuous(limits = requested_sens_ylim, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.margin = margin(34, 0, 48, 16))
fig_requested_sens_forest <- ggplot(requested_sens_fig_df, aes(y = row, x = HR)) +
  geom_hline(aes(yintercept = row), color = "#E6E6E6", linewidth = 0.55) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.8, color = "#555555") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                width = 0.12, linewidth = 0.95, color = COL_BLACK) +
  geom_point(size = 4.4, color = COL_BLACK) +
  scale_y_continuous(limits = requested_sens_ylim, breaks = NULL,
                     labels = NULL, expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.74, 1.04), breaks = seq(0.75, 1.00, 0.05)) +
  labs(title = "Sensitivity analyses",
       x = "Hazard ratio (ultra-early vs early)", y = NULL) +
  theme_icm(18) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 21, margin = margin(t = 10)),
        plot.title = element_text(size = 25, face = "bold", hjust = 0),
        plot.margin = margin(10, 0, 18, 0))
fig_requested_sens_values <- ggplot(requested_sens_fig_df, aes(x = 0, y = row, label = text)) +
  geom_text(hjust = 0, size = 4.2, family = FONT_FAMILY, fontface = "bold",
            color = COL_BLACK) +
  scale_y_continuous(limits = requested_sens_ylim, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.margin = margin(34, 18, 48, 6))
fig_requested_sens <- fig_requested_sens_labels + fig_requested_sens_forest + fig_requested_sens_values +
  plot_layout(widths = c(0.42, 0.42, 0.16))

vp_dose_ext_table <- pt_ne %>%
  filter(!is.na(tvp_rel)) %>%
  mutate(`Observed VP timing` = make_timing_label(obs_vp_group)) %>%
  filter(!is.na(`Observed VP timing`), `Observed VP timing` != "Non-VP use") %>%
  group_by(`Observed VP timing`) %>%
  summarise(
    n = n(),
    `Initial dose, U/min` = fmt_med_iqr(vp_dose_init, 3),
    `Maximum dose within 6 h, U/min` = fmt_med_iqr(vp_dose_max_6h, 3),
    `Maximum dose within 24 h, U/min` = fmt_med_iqr(vp_dose_max_24h, 3),
    `Cumulative dose within 6 h, U` = fmt_med_iqr(vp_cum_6h, 2),
    `Cumulative dose within 24 h, U` = fmt_med_iqr(vp_cum_24h, 2),
    `Initial dose <0.03 U/min` = fmt_count_pct(coalesce(vp_init_lt003, FALSE)),
    `Initial dose =0.03 U/min` = fmt_count_pct(coalesce(vp_init_eq003, FALSE)),
    `Initial dose >0.03 U/min` = fmt_count_pct(coalesce(vp_init_gt003, FALSE)),
    `Dose data missing` = fmt_count_pct(coalesce(vp_dose_missing, FALSE)),
    .groups = "drop"
  ) %>%
  mutate(n = fmt_n(n))
ft_vp_dose_ext <- make_ft(
  vp_dose_ext_table,
  "Supplementary Table. Vasopressin dose by observed timing",
  8.0
) %>%
  add_ft_note("Values are median (IQR) or n (%), restricted to patients with observed VP initiation.", 7.0)

transpose_group_table <- function(tbl, group_col, metric_col = "Characteristic") {
  groups <- as.character(tbl[[group_col]])
  metric_names <- setdiff(names(tbl), group_col)
  map_dfr(metric_names, function(metric) {
    vals <- as.list(tbl[[metric]])
    names(vals) <- groups
    bind_cols(tibble(!!metric_col := metric), as_tibble(vals, .name_repair = "minimal"))
  })
}
vp_dose_bundle_table <- transpose_group_table(vp_dose_ext_table, "Observed VP timing")
ft_vp_dose_bundle <- make_ft(
  vp_dose_bundle_table,
  "Supplementary Table. Vasopressin dose by observed timing",
  8.4
) %>%
  add_ft_note("Values are median (IQR) or n (%), restricted to patients with observed VP initiation.", 7.1)

vp_obs_for_traj <- pt_ne %>%
  filter(!is.na(tvp_rel), tvp_rel <= EARLY_END) %>%
  transmute(uid, t0_hr, tvp_rel,
            timing = factor(if_else(tvp_rel <= 3, "Ultra-early (0-3 h)", "Early (>3-6 h)"),
                            levels = c("Ultra-early (0-3 h)", "Early (>3-6 h)")))

assigned_adherent_weights <- ccw_tbl_main %>%
  mutate(
    w = pmax(P_CLIP, pmin(weight, caps_main$cap99)),
    timing = case_when(
      arm == "ultra" & !is.na(tvp_rel) & tvp_rel <= 3 ~ "Ultra-early (0-3 h)",
      arm == "early" & !is.na(tvp_rel) & tvp_rel > 3 & tvp_rel <= 6 ~ "Early (>3-6 h)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(timing)) %>%
  transmute(uid, t0_hr, tvp_rel, timing = factor(timing, levels = c("Ultra-early (0-3 h)", "Early (>3-6 h)")), w)

weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
  ok <- !is.na(x) & !is.na(w) & is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

fmt_wmed_iqr <- function(x, w, digits = 1) {
  q <- weighted_quantile(x, w, c(0.25, 0.5, 0.75))
  if (any(is.na(q))) return("")
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"), q[2], q[1], q[3])
}

fmt_event_wpct <- function(event, w) {
  ok <- !is.na(event) & !is.na(w) & is.finite(w) & w > 0
  if (!any(ok)) return("")
  sprintf("%s; %.1f%%", fmt_n(sum(event[ok], na.rm = TRUE)),
          100 * weighted.mean(as.numeric(event[ok]), w[ok], na.rm = TRUE))
}

lactate_trajectory_patient <- ts_imputed %>%
  inner_join(assigned_adherent_weights, by = "uid") %>%
  mutate(t_rel_vp = time_hr - (t0_hr + tvp_rel),
         hour_after_vp = as.integer(round(t_rel_vp))) %>%
  filter(hour_after_vp %in% c(0L, 3L, 6L, 12L, 24L)) %>%
  group_by(uid, timing, w, hour_after_vp) %>%
  summarise(lactate = first_nonmissing(lact), .groups = "drop")

lactate_trajectory_summary <- lactate_trajectory_patient %>%
  group_by(timing, hour_after_vp) %>%
  summarise(
    n = sum(!is.na(lactate)),
    q = list(weighted_quantile(lactate, w, c(0.25, 0.5, 0.75))),
    .groups = "drop"
  ) %>%
  mutate(q1 = map_dbl(q, 1), median = map_dbl(q, 2), q3 = map_dbl(q, 3)) %>%
  select(-q) %>%
  mutate(`Median (IQR), mmol/L` = sprintf("%.1f (%.1f, %.1f)", median, q1, q3))

lactate_trajectory_table <- lactate_trajectory_summary %>%
  transmute(`Assigned strategy` = timing,
            `Hours after VP initiation` = hour_after_vp,
            `Patients with lactate, n` = fmt_n(n),
            `Weighted lactate, mmol/L` = `Median (IQR), mmol/L`)
ft_lactate_traj <- make_ft(
  lactate_trajectory_table,
  "Supplementary Table. IPCW-weighted lactate trajectory after vasopressin initiation",
  8.2
)

fig_lactate_traj <- ggplot(lactate_trajectory_summary,
                           aes(x = hour_after_vp, y = median, linetype = timing, group = timing)) +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.55, linewidth = 0.85, color = COL_BLACK) +
  geom_line(linewidth = 1.15, color = COL_BLACK) +
  geom_point(size = 3.8, color = COL_BLACK) +
  geom_text(aes(label = fmt_n(n)), vjust = -1.15, size = 4.0,
            family = FONT_FAMILY, fontface = "bold", color = COL_BLACK,
            show.legend = FALSE) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_x_continuous(breaks = c(0, 3, 6, 12, 24)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  labs(x = "Hours after VP initiation",
       y = "Weighted lactate, mmol/L, median (IQR)",
       linetype = NULL) +
  theme_icm(18) +
  theme(legend.position = c(0.73, 0.86),
        legend.background = element_rect(fill = alpha("white", 0.90), color = NA),
        axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)))

post_vp_agent_flags <- ts_imputed %>%
  inner_join(pt_ne %>% filter(!is.na(tvp_rel)) %>% select(uid, t0_hr, tvp_rel),
             by = "uid") %>%
  mutate(t_rel_vp = time_hr - (t0_hr + tvp_rel)) %>%
  group_by(uid) %>%
  summarise(
    dopa_pre_or_at = any(replace_na(dopa_rate, 0) > 0 & t_rel_vp <= 0 & t_rel_vp >= -24),
    dobu_pre_or_at = any(replace_na(dobu_rate, 0) > 0 & t_rel_vp <= 0 & t_rel_vp >= -24),
    epi_pre_or_at = any(replace_na(epi_rate, 0) > 0 & t_rel_vp <= 0 & t_rel_vp >= -24),
    new_dopa_3h = !dopa_pre_or_at & any(replace_na(dopa_rate, 0) > 0 & t_rel_vp > 0 & t_rel_vp <= 3),
    new_dopa_6h = !dopa_pre_or_at & any(replace_na(dopa_rate, 0) > 0 & t_rel_vp > 0 & t_rel_vp <= 6),
    new_dopa_24h = !dopa_pre_or_at & any(replace_na(dopa_rate, 0) > 0 & t_rel_vp > 0 & t_rel_vp <= 24),
    new_dobu_3h_derived = !dobu_pre_or_at & any(replace_na(dobu_rate, 0) > 0 & t_rel_vp > 0 & t_rel_vp <= 3),
    new_dobu_6h_derived = !dobu_pre_or_at & any(replace_na(dobu_rate, 0) > 0 & t_rel_vp > 0 & t_rel_vp <= 6),
    new_dobu_24h_derived = !dobu_pre_or_at & any(replace_na(dobu_rate, 0) > 0 & t_rel_vp > 0 & t_rel_vp <= 24),
    new_epi_3h_derived = !epi_pre_or_at & any(replace_na(epi_rate, 0) > 0 & t_rel_vp > 0 & t_rel_vp <= 3),
    new_epi_6h_derived = !epi_pre_or_at & any(replace_na(epi_rate, 0) > 0 & t_rel_vp > 0 & t_rel_vp <= 6),
    new_epi_24h_derived = !epi_pre_or_at & any(replace_na(epi_rate, 0) > 0 & t_rel_vp > 0 & t_rel_vp <= 24),
    .groups = "drop"
  )

post_vp_events_table <- pt_ne %>%
  inner_join(assigned_adherent_weights %>% select(uid, timing, w), by = "uid") %>%
  left_join(post_vp_agent_flags, by = "uid") %>%
  mutate(`Assigned strategy` = timing) %>%
  group_by(`Assigned strategy`) %>%
  summarise(
    `IPCW analytic clones, n` = n(),
    `VP interruption within 3 h, n; weighted %` = fmt_event_wpct(coalesce(interrupt_3h, FALSE), w),
    `Permanent VP discontinuation within 3 h, n; weighted %` = fmt_event_wpct(coalesce(permanent_disc_3h, FALSE), w),
    `New dopamine within 3 h, n; weighted %` = fmt_event_wpct(coalesce(new_dopa_3h, FALSE), w),
    `New dopamine within 6 h, n; weighted %` = fmt_event_wpct(coalesce(new_dopa_6h, FALSE), w),
    `New dopamine within 24 h, n; weighted %` = fmt_event_wpct(coalesce(new_dopa_24h, FALSE), w),
    `New dobutamine within 3 h, n; weighted %` = fmt_event_wpct(coalesce(new_dobu_3h_derived, FALSE), w),
    `New dobutamine within 6 h, n; weighted %` = fmt_event_wpct(coalesce(new_dobu_6h_derived, FALSE), w),
    `New dobutamine within 24 h, n; weighted %` = fmt_event_wpct(coalesce(new_dobu_24h_derived, FALSE), w),
    `New epinephrine within 3 h, n; weighted %` = fmt_event_wpct(coalesce(new_epi_3h_derived, FALSE), w),
    `New epinephrine within 6 h, n; weighted %` = fmt_event_wpct(coalesce(new_epi_6h_derived, FALSE), w),
    `New epinephrine within 24 h, n; weighted %` = fmt_event_wpct(coalesce(new_epi_24h_derived, FALSE), w),
    .groups = "drop"
  ) %>%
  mutate(`IPCW analytic clones, n` = fmt_n(`IPCW analytic clones, n`))
ft_post_vp_events <- make_ft(
  post_vp_events_table,
  "Supplementary Table. IPCW-weighted cardiovascular support signals after vasopressin",
  7.1
) %>%
  add_ft_note(
    "Values are unweighted event counts and IPCW-weighted percentages among clones adherent to the assigned strategy. New dopamine, dobutamine, and epinephrine indicate no use during the 24 h before or at VP initiation and new use after VP initiation.",
    6.6
  )
post_vp_events_bundle_table <- transpose_group_table(post_vp_events_table, "Assigned strategy")
ft_post_vp_events_bundle <- make_ft(
  post_vp_events_bundle_table,
  "Supplementary Table. IPCW-weighted post-vasopressin interruption and new vasoactive/inotropic support",
  8.2
) %>%
  add_ft_note(
    "Values are unweighted event counts and IPCW-weighted percentages among clones adherent to the assigned strategy. New dopamine, dobutamine, and epinephrine indicate no use during the 24 h before or at VP initiation and new use after VP initiation.",
    7.0
  )

safety_desc_table <- pt_ne %>%
  mutate(`Observed VP timing` = make_timing_label(obs_vp_group)) %>%
  group_by(`Observed VP timing`) %>%
  summarise(
    n = n(),
    DVT = fmt_count_pct(coalesce(sf_dvt, FALSE)),
    PE = fmt_count_pct(coalesce(sf_pe, FALSE)),
    `Mesenteric ischemia` = fmt_count_pct(coalesce(sf_mesenteric, FALSE)),
    `Digital/peripheral ischemia` = fmt_count_pct(coalesce(sf_digital, FALSE)),
    .groups = "drop"
  ) %>%
  mutate(n = fmt_n(n))
ft_safety_desc <- make_ft(
  safety_desc_table,
  "Supplementary Table. Exploratory safety outcomes by observed timing",
  8.5
) %>%
  add_ft_note("Values are n (%). These diagnosis-based safety outcomes are exploratory.", 7.2)
safety_desc_bundle_table <- transpose_group_table(safety_desc_table, "Observed VP timing")
ft_safety_desc_bundle <- make_ft(
  safety_desc_bundle_table,
  "Supplementary Table. Exploratory safety outcomes by observed timing",
  8.5
) %>%
  add_ft_note("Values are n (%). These diagnosis-based safety outcomes are exploratory.", 7.2)

fit_weighted_poisson_binary <- function(ccw_tbl, pt, outcome, label) {
  dat <- ccw_tbl %>%
    select(uid, arm, src, weight) %>%
    left_join(pt %>% select(uid, all_of(outcome)), by = "uid") %>%
    mutate(
      .outcome = coalesce(as.logical(.data[[outcome]]), FALSE),
      w = pmax(P_CLIP, pmin(weight, caps_main$cap99)),
      arm = factor(arm, levels = c("early", "ultra")),
      src = factor(src)
    )
  event_by_uid <- dat %>%
    group_by(uid) %>%
    summarise(.event = any(.outcome, na.rm = TRUE), .groups = "drop")
  if (sum(dat$.outcome, na.rm = TRUE) < 2) {
    return(tibble(label = label, RR = NA_real_, lower = NA_real_, upper = NA_real_,
                  n_patients = n_distinct(dat$uid), n_event = sum(event_by_uid$.event, na.rm = TRUE)))
  }
  fit <- suppressWarnings(glm(.outcome ~ arm + src, family = poisson(link = "log"),
                              data = dat, weights = w))
  vc <- sandwich::vcovCL(fit, cluster = dat$uid)
  est <- coef(fit)["armultra"]
  se <- sqrt(vc["armultra", "armultra"])
  tibble(
    label = label,
    RR = exp(est),
    lower = exp(est - 1.96 * se),
    upper = exp(est + 1.96 * se),
    n_patients = n_distinct(dat$uid),
    n_event = sum(event_by_uid$.event, na.rm = TRUE)
  )
}

safety_model_results <- bind_rows(
  fit_weighted_poisson_binary(ccw_tbl_main, pt_ne, "sf_dvt", "DVT"),
  fit_weighted_poisson_binary(ccw_tbl_main, pt_ne, "sf_pe", "PE"),
  fit_weighted_poisson_binary(ccw_tbl_main, pt_ne, "sf_mesenteric", "Mesenteric ischemia"),
  fit_weighted_poisson_binary(ccw_tbl_main, pt_ne, "sf_digital", "Digital/peripheral ischemia")
) %>%
  mutate(`RR (95% CI)` = fmt_hr(RR, lower, upper),
         `Model patients, n` = fmt_n(n_patients),
         `Events, n` = fmt_n(n_event))
safety_model_table <- safety_model_results %>%
  transmute(Outcome = label, `Model patients, n`, `Events, n`, `RR (95% CI)`)
ft_safety_model <- make_ft(
  safety_model_table,
  "Supplementary Table. Exploratory safety outcome models",
  8.8
) %>%
  add_ft_note(
    "Risk ratios compare the ultra-early strategy with the early strategy using IPCW-weighted Poisson models with robust clustered standard errors. Safety outcomes are diagnosis based and exploratory.",
    7.1
  )

survivor_bundle_table <- transpose_group_table(survivor_table, "28-d status")
ft_survivor_bundle <- make_ft(
  survivor_bundle_table,
  "Supplementary Table. Survivors versus non-survivors at 28 days",
  8.5
) %>%
  add_ft_note("Values are median (IQR) or n (%).", 7.2)

make_rr_forest <- function(df, title = NULL) {
  plot_df <- df %>%
    mutate(row = rev(row_number()), text = fmt_hr(RR, lower, upper))
  x_low <- min(0.95, max(0.10, min(plot_df$lower, na.rm = TRUE) * 0.80))
  x_high <- max(1.05, min(10, max(plot_df$upper, na.rm = TRUE) * 1.25))
  breaks <- c(0.25, 0.5, 1, 2, 4, 8)
  ggplot(plot_df, aes(y = row, x = RR)) +
    geom_vline(xintercept = 1, linetype = "22", linewidth = 0.8, color = "#555555") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                  width = 0.12, linewidth = 0.95, color = COL_BLACK) +
    geom_point(size = 4.6, color = COL_BLACK) +
    geom_text(aes(x = x_high / 1.05, label = text),
              hjust = 1, size = 4.8, family = FONT_FAMILY, fontface = "bold") +
    scale_y_continuous(breaks = plot_df$row, labels = plot_df$label,
                       expand = expansion(add = 0.6)) +
    scale_x_log10(limits = c(x_low, x_high), breaks = breaks,
                  labels = function(x) sprintf("%g", x)) +
    labs(title = title, x = "Risk ratio (ultra-early vs early)", y = NULL) +
    theme_icm(18) +
    theme(axis.text.y = element_text(size = 16, face = "bold", hjust = 1),
          axis.title.x = element_text(size = 21, margin = margin(t = 10)),
          plot.title = element_text(size = 25, face = "bold"))
}
safety_model_plot_df <- safety_model_results %>%
  transmute(label, RR, lower, upper, text = fmt_hr(RR, lower, upper)) %>%
  mutate(row = rev(row_number()))
safety_model_ylim <- c(0.5, nrow(safety_model_plot_df) + 0.5)
safety_x_low <- min(0.95, max(0.10, min(safety_model_plot_df$lower, na.rm = TRUE) * 0.80))
safety_x_high <- max(1.05, min(10, max(safety_model_plot_df$upper, na.rm = TRUE) * 1.25))
fig_safety_labels <- ggplot(safety_model_plot_df, aes(y = row, label = label)) +
  geom_text(aes(x = 0), hjust = 0, size = 4.65, family = FONT_FAMILY,
            fontface = "bold", color = COL_BLACK) +
  scale_y_continuous(limits = safety_model_ylim, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.margin = margin(34, 0, 48, 16))
fig_safety_forest <- ggplot(safety_model_plot_df, aes(y = row, x = RR)) +
  geom_hline(aes(yintercept = row), color = "#E6E6E6", linewidth = 0.55) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.8, color = "#555555") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                width = 0.12, linewidth = 0.95, color = COL_BLACK) +
  geom_point(size = 4.4, color = COL_BLACK) +
  scale_y_continuous(limits = safety_model_ylim, breaks = NULL,
                     labels = NULL, expand = c(0, 0)) +
  scale_x_log10(limits = c(safety_x_low, safety_x_high),
                breaks = c(0.25, 0.5, 1, 2, 4, 8),
                labels = function(x) sprintf("%g", x)) +
  labs(title = "Exploratory safety outcomes",
       x = "Risk ratio (ultra-early vs early)", y = NULL) +
  theme_icm(18) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 21, margin = margin(t = 10)),
        plot.title = element_text(size = 25, face = "bold", hjust = 0),
        plot.margin = margin(10, 0, 18, 0))
fig_safety_values <- ggplot(safety_model_plot_df, aes(x = 0, y = row, label = text)) +
  geom_text(hjust = 0, size = 4.2, family = FONT_FAMILY, fontface = "bold",
            color = COL_BLACK) +
  scale_y_continuous(limits = safety_model_ylim, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.margin = margin(34, 18, 48, 6))
fig_safety_model <- fig_safety_labels + fig_safety_forest + fig_safety_values +
  plot_layout(widths = c(0.34, 0.48, 0.18))

OR_REVISED_DIR <- file.path(ROOT_DIR, "Outcome_revise後", "OR7_16_revised_assets")
dir_create(OR_REVISED_DIR)

assigned_strategy_weights <- assigned_adherent_weights %>%
  group_by(uid, timing) %>%
  summarise(t0_hr = first(t0_hr), tvp_rel = first(tvp_rel), w = last(w), .groups = "drop")

save_or_figure <- function(plot_obj, stem, width = 11, height = 7.5) {
  ggsave(file.path(OR_REVISED_DIR, paste0(stem, ".pdf")),
         plot_obj, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(OR_REVISED_DIR, paste0(stem, ".png")),
         plot_obj, width = width, height = height, dpi = 220, bg = "white")
}

make_clean_forest <- function(df, title, xlim = c(0.68, 1.16), breaks = seq(0.7, 1.1, 0.1)) {
  plot_df <- df %>%
    mutate(row = rev(row_number()), text = fmt_hr(HR, lower, upper))
  ggplot(plot_df, aes(y = row, x = HR)) +
    geom_vline(xintercept = 1, linetype = "22", linewidth = 0.8, color = "#555555") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                  width = 0.12, linewidth = 0.95, color = COL_BLACK) +
    geom_point(size = 4.4, color = COL_BLACK) +
    geom_text(aes(x = xlim[2] - 0.01, label = text),
              hjust = 1, size = 4.6, family = FONT_FAMILY, fontface = "bold") +
    scale_y_continuous(breaks = plot_df$row, labels = plot_df$label,
                       expand = expansion(add = 0.6)) +
    scale_x_continuous(limits = xlim, breaks = breaks) +
    labs(title = title, x = "Hazard ratio (ultra-early vs early)", y = NULL) +
    theme_icm(18) +
    theme(axis.text.y = element_text(size = 15.5, face = "bold", hjust = 1),
          axis.title.x = element_text(size = 21, margin = margin(t = 10)),
          plot.title = element_text(size = 25, face = "bold", hjust = 0))
}

fig_or7_vp_timing <- vp_timing_bins %>%
  mutate(strategy = if_else(as.integer(timing_bin) <= 3,
                            "Ultra-early (0-3 h)", "Early (>3-6 h)")) %>%
  ggplot(aes(x = timing_bin, y = pct * 100, fill = strategy)) +
  geom_col(color = COL_BLACK, linewidth = 0.85, width = 0.72) +
  geom_text(aes(label = paste0(fmt_n(n), "\n", sprintf("%.1f%%", pct * 100))),
            vjust = -0.18, size = 4.9, family = FONT_FAMILY, fontface = "bold") +
  geom_vline(xintercept = 3.5, linetype = "22", linewidth = 0.75, color = "#555555") +
  scale_fill_manual(values = c("Ultra-early (0-3 h)" = "white",
                               "Early (>3-6 h)" = "#D9D9D9")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Hours from NE time zero to first VP initiation",
       y = "Patients with VP initiation within 6 h (%)",
       fill = NULL) +
  theme_icm(18) +
  theme(legend.position = c(0.74, 0.84),
        legend.background = element_rect(fill = alpha("white", 0.90), color = NA),
        axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)))

or8_weight_plot_df <- ccw_tbl_main %>%
  mutate(
    `Assigned strategy` = factor(arm, levels = c("ultra", "early"),
                                 labels = c("Ultra-early (0-3 h)", "Early (>3-6 h)")),
    `Truncated IPCW` = pmax(P_CLIP, pmin(weight, caps_main$cap99))
  )
or8_weight_label <- or8_weight_plot_df %>%
  group_by(`Assigned strategy`) %>%
  summarise(
    med = median(`Truncated IPCW`, na.rm = TRUE),
    q1 = quantile(`Truncated IPCW`, 0.25, na.rm = TRUE),
    q3 = quantile(`Truncated IPCW`, 0.75, na.rm = TRUE),
    max_w = max(`Truncated IPCW`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("Median %.2f (IQR %.2f-%.2f)", med, q1, q3))
fig_or8_weights <- ggplot(or8_weight_plot_df, aes(x = `Truncated IPCW`)) +
  geom_histogram(bins = 38, fill = "white", color = COL_BLACK, linewidth = 0.65) +
  geom_vline(data = or8_weight_label, aes(xintercept = med),
             linetype = "22", linewidth = 0.8, color = "#555555") +
  geom_text(data = or8_weight_label, aes(x = max_w, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 1.02, vjust = 1.5, size = 4.6,
            family = FONT_FAMILY, fontface = "bold") +
  facet_wrap(~`Assigned strategy`, ncol = 1, scales = "free_y") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20),
                labels = c("0.5", "1", "2", "5", "10", "20")) +
  labs(x = "Stabilized IPCW after 99th-percentile truncation, log scale",
       y = "Clone intervals, n") +
  theme_icm(18) +
  theme(strip.text = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)))

fig_or10_database <- make_clean_forest(
  database_effects %>% select(label, HR, lower, upper),
  "Database-specific primary outcome estimates",
  xlim = c(0.68, 1.16)
)

or11_plot_df <- or11_sensitivity_results %>%
  mutate(
    label = factor(label, levels = c(
      "Primary: NE time zero",
      "Primary: 95th percentile truncation",
      "Database: MIMIC-IV",
      "Database: eICU-CRD",
      "Sensitivity: NEE time zero",
      "Sensitivity: Nad-only time zero",
      "Sensitivity: excluding dobutamine at time zero"
    ))
  ) %>%
  arrange(label) %>%
  mutate(label = as.character(label)) %>%
  select(label, HR, lower, upper)
fig_or11_sensitivity <- make_clean_forest(
  or11_plot_df,
  "Primary outcome sensitivity analyses",
  xlim = c(0.68, 1.16)
)

nee_traj_patient <- ts_imputed %>%
  inner_join(assigned_strategy_weights, by = "uid") %>%
  mutate(
    hour_after_t0 = as.integer(round(time_hr - t0_hr)),
    nee = coalesce(norepi_equiv, norepi_rate, 0)
  ) %>%
  filter(hour_after_t0 >= 0, hour_after_t0 <= 72) %>%
  group_by(uid, timing, w, hour_after_t0) %>%
  summarise(nee = first_nonmissing(nee), .groups = "drop")

nee_traj_summary <- nee_traj_patient %>%
  group_by(timing, hour_after_t0) %>%
  summarise(
    n = sum(!is.na(nee)),
    q = list(weighted_quantile(nee, w, c(0.25, 0.5, 0.75))),
    .groups = "drop"
  ) %>%
  mutate(q1 = map_dbl(q, 1), median = map_dbl(q, 2), q3 = map_dbl(q, 3)) %>%
  select(-q)
nee_traj_table <- nee_traj_summary %>%
  filter(hour_after_t0 %in% c(0, 3, 6, 12, 24, 48, 72)) %>%
  transmute(
    `Assigned strategy` = timing,
    `Hours after time zero` = hour_after_t0,
    `Patients with data, n` = fmt_n(n),
    `Weighted NEE dose, mcg/kg/min` = sprintf("%.2f (%.2f, %.2f)", median, q1, q3)
  )
fig_or13_nee_traj <- ggplot(nee_traj_summary,
                            aes(x = hour_after_t0, y = median, linetype = timing, group = timing)) +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = timing),
              alpha = 0.16, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2, color = COL_BLACK) +
  geom_point(data = nee_traj_summary %>% filter(hour_after_t0 %% 12 == 0),
             size = 3.0, color = COL_BLACK) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_fill_manual(values = c("Ultra-early (0-3 h)" = "#A0A0A0",
                               "Early (>3-6 h)" = "#D0D0D0")) +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Hours after NE time zero",
       y = "Weighted NEE dose, mcg/kg/min, median (IQR)",
       linetype = NULL) +
  theme_icm(18) +
  theme(legend.position = c(0.73, 0.86),
        legend.background = element_rect(fill = alpha("white", 0.90), color = NA),
        axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)))

weighted_sd <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (!length(x)) return(NA_real_)
  mu <- weighted.mean(x, w)
  sqrt(sum(w * (x - mu)^2) / sum(w))
}
fmt_wmed_iqr_bracket <- function(x, w, digits = 1) {
  q <- weighted_quantile(x, w, c(0.25, 0.5, 0.75))
  if (any(is.na(q))) return("")
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), q[2], q[1], q[3])
}
fmt_wmean_sd <- function(x, w, digits = 2) {
  mu <- weighted.mean(x, w, na.rm = TRUE)
  sd <- weighted_sd(x, w)
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mu, sd)
}
weighted_diff_ci <- function(df, var) {
  dat <- df %>%
    mutate(strategy = factor(timing, levels = c("Early (>3-6 h)", "Ultra-early (0-3 h)"))) %>%
    filter(!is.na(.data[[var]]), !is.na(w), is.finite(w), w > 0)
  fit <- lm(as.formula(paste(var, "~ strategy")), data = dat, weights = w)
  vc <- sandwich::vcovHC(fit, type = "HC0")
  est <- unname(coef(fit)[2])
  se <- sqrt(vc[2, 2])
  sprintf("%.2f (%.2f to %.2f)", est, est - 1.96 * se, est + 1.96 * se)
}

post_t0_exposure <- ts_imputed %>%
  inner_join(assigned_strategy_weights, by = "uid") %>%
  mutate(
    t_rel = time_hr - t0_hr,
    nee = coalesce(norepi_equiv, norepi_rate, 0),
    any_vasopressor = coalesce(norepi_equiv, 0) > 0 |
      coalesce(norepi_rate, 0) > 0 | coalesce(epi_rate, 0) > 0 |
      coalesce(dopa_rate, 0) > 0 | coalesce(phn_rate, 0) > 0 |
      coalesce(adh_rate, 0) > 0 | coalesce(vaso_ind, FALSE)
  ) %>%
  filter(t_rel >= 0, t_rel <= DAY28_HR)

exposure_patient <- post_t0_exposure %>%
  group_by(uid, timing, w) %>%
  summarise(
    nee_auc_72h = sum(nee[t_rel <= 72], na.rm = TRUE),
    nee_max_72h = if (all(is.na(nee[t_rel <= 72]))) NA_real_ else max(nee[t_rel <= 72], na.rm = TRUE),
    vp_hrs_7d = sum(any_vasopressor[t_rel <= 7 * 24], na.rm = TRUE),
    vp_hrs_14d = sum(any_vasopressor[t_rel <= 14 * 24], na.rm = TRUE),
    vp_hrs_28d = sum(any_vasopressor[t_rel <= DAY28_HR], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pt_ne %>% select(uid, death_time, follow_end_rel), by = "uid") %>%
  mutate(
    follow_7 = pmin(coalesce(follow_end_rel, 7 * 24), 7 * 24),
    follow_14 = pmin(coalesce(follow_end_rel, 14 * 24), 14 * 24),
    follow_28 = pmin(coalesce(follow_end_rel, DAY28_HR), DAY28_HR),
    died_7d = is.finite(death_time) & death_time <= 7 * 24,
    died_14d = is.finite(death_time) & death_time <= 14 * 24,
    died_28d = is.finite(death_time) & death_time <= DAY28_HR,
    vfd_7 = if_else(died_7d, 0, pmax(follow_7 - vp_hrs_7d, 0) / 24),
    vfd_14 = if_else(died_14d, 0, pmax(follow_14 - vp_hrs_14d, 0) / 24),
    vfd_28 = if_else(died_28d, 0, pmax(follow_28 - vp_hrs_28d, 0) / 24)
  )

or14_table <- tibble(
  Variable = c(
    "NEE AUC within 72 h, mcg/kg/min-hours",
    "Maximum NEE dose within 72 h, mcg/kg/min",
    "Vasopressor-free days to day 7, mean (SD)",
    "Vasopressor-free days to day 14, mean (SD)",
    "Vasopressor-free days to day 28, mean (SD)"
  ),
  `Ultra-early (0-3 h)` = c(
    fmt_wmed_iqr_bracket(exposure_patient$nee_auc_72h[exposure_patient$timing == "Ultra-early (0-3 h)"],
                         exposure_patient$w[exposure_patient$timing == "Ultra-early (0-3 h)"], 1),
    fmt_wmed_iqr_bracket(exposure_patient$nee_max_72h[exposure_patient$timing == "Ultra-early (0-3 h)"],
                         exposure_patient$w[exposure_patient$timing == "Ultra-early (0-3 h)"], 2),
    fmt_wmean_sd(exposure_patient$vfd_7[exposure_patient$timing == "Ultra-early (0-3 h)"],
                 exposure_patient$w[exposure_patient$timing == "Ultra-early (0-3 h)"], 2),
    fmt_wmean_sd(exposure_patient$vfd_14[exposure_patient$timing == "Ultra-early (0-3 h)"],
                 exposure_patient$w[exposure_patient$timing == "Ultra-early (0-3 h)"], 2),
    fmt_wmean_sd(exposure_patient$vfd_28[exposure_patient$timing == "Ultra-early (0-3 h)"],
                 exposure_patient$w[exposure_patient$timing == "Ultra-early (0-3 h)"], 2)
  ),
  `Early (>3-6 h)` = c(
    fmt_wmed_iqr_bracket(exposure_patient$nee_auc_72h[exposure_patient$timing == "Early (>3-6 h)"],
                         exposure_patient$w[exposure_patient$timing == "Early (>3-6 h)"], 1),
    fmt_wmed_iqr_bracket(exposure_patient$nee_max_72h[exposure_patient$timing == "Early (>3-6 h)"],
                         exposure_patient$w[exposure_patient$timing == "Early (>3-6 h)"], 2),
    fmt_wmean_sd(exposure_patient$vfd_7[exposure_patient$timing == "Early (>3-6 h)"],
                 exposure_patient$w[exposure_patient$timing == "Early (>3-6 h)"], 2),
    fmt_wmean_sd(exposure_patient$vfd_14[exposure_patient$timing == "Early (>3-6 h)"],
                 exposure_patient$w[exposure_patient$timing == "Early (>3-6 h)"], 2),
    fmt_wmean_sd(exposure_patient$vfd_28[exposure_patient$timing == "Early (>3-6 h)"],
                 exposure_patient$w[exposure_patient$timing == "Early (>3-6 h)"], 2)
  ),
  `Difference (95% CI)` = c(
    "",
    "",
    weighted_diff_ci(exposure_patient, "vfd_7"),
    weighted_diff_ci(exposure_patient, "vfd_14"),
    weighted_diff_ci(exposure_patient, "vfd_28")
  )
)
ft_or14 <- make_ft(or14_table, "Online Resource 14. Weighted catecholamine exposure and vasopressor-free days", 8.5) %>%
  add_ft_note(
    "AUC is calculated as the sum of hourly norepinephrine-equivalent dose from time zero through 72 h. Vasopressor-free days count days alive and free from any vasopressor; patients who died before each horizon were assigned 0 vasopressor-free days.",
    7.0
  )

or15_dist_df <- pt_ne %>%
  filter(obs_vp_group %in% c("0-3 h", ">3-6 h")) %>%
  transmute(
    `Observed VP timing` = factor(make_timing_label(obs_vp_group),
                                  levels = c("Ultra-early (0-3 h)", "Early (>3-6 h)")),
    `Creatinine at time zero, mg/dL` = crea0,
    `Lactate at time zero, mmol/L` = lact0,
    `MAP at time zero, mmHg` = map0,
    `Norepinephrine dose at time zero, mcg/kg/min` = ne_rate_0
  ) %>%
  pivot_longer(-`Observed VP timing`, names_to = "Variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(Variable) %>%
  filter(value <= quantile(value, 0.99, na.rm = TRUE)) %>%
  ungroup()
fig_or15_baseline_dist <- ggplot(or15_dist_df,
                                 aes(x = value, linetype = `Observed VP timing`,
                                     fill = `Observed VP timing`)) +
  geom_density(alpha = 0.15, color = COL_BLACK, linewidth = 1.0) +
  facet_wrap(~Variable, scales = "free", ncol = 2) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_fill_manual(values = c("Ultra-early (0-3 h)" = "#A0A0A0",
                               "Early (>3-6 h)" = "#D0D0D0")) +
  labs(x = NULL, y = "Density", linetype = NULL, fill = NULL) +
  theme_icm(17) +
  theme(legend.position = "top",
        strip.text = element_text(size = 15.5, face = "bold"),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)))

hemo_traj_summary <- ts_imputed %>%
  inner_join(pt_ne %>%
               filter(obs_vp_group %in% c("0-3 h", ">3-6 h")) %>%
               transmute(uid, t0_hr,
                         `Observed VP timing` = factor(make_timing_label(obs_vp_group),
                                                       levels = c("Ultra-early (0-3 h)", "Early (>3-6 h)"))),
             by = "uid") %>%
  mutate(hour = as.integer(round(time_hr - t0_hr))) %>%
  filter(hour >= -6, hour <= 6) %>%
  transmute(
    uid, `Observed VP timing`, hour,
    `MAP, mmHg` = map,
    `Lactate, mmol/L` = lact,
    `NEE dose, mcg/kg/min` = coalesce(norepi_equiv, norepi_rate, 0)
  ) %>%
  pivot_longer(c(`MAP, mmHg`, `Lactate, mmol/L`, `NEE dose, mcg/kg/min`),
               names_to = "Variable", values_to = "value") %>%
  group_by(`Observed VP timing`, hour, Variable) %>%
  summarise(
    n = sum(!is.na(value)),
    median = median(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
fig_or16_hemo_traj <- ggplot(hemo_traj_summary,
                             aes(x = hour, y = median, linetype = `Observed VP timing`,
                                 group = `Observed VP timing`)) +
  geom_line(linewidth = 1.1, color = COL_BLACK) +
  geom_point(size = 2.8, color = COL_BLACK) +
  geom_vline(xintercept = 0, linetype = "11", linewidth = 0.65, color = "#555555") +
  facet_wrap(~Variable, scales = "free_y", ncol = 1) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_x_continuous(breaks = seq(-6, 6, 3)) +
  labs(x = "Hours relative to NE time zero", y = "Median value", linetype = NULL) +
  theme_icm(17) +
  theme(legend.position = "top",
        strip.text = element_text(size = 15.5, face = "bold"),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)))

ft_or10_db <- make_ft(
  database_effects %>%
    transmute(Database = label,
              `Clones, n` = fmt_n(n_clone),
              `Events, n` = fmt_n(n_event),
              `HR (95% CI)` = `HR (95% CI)`),
  "Online Resource 10. Database-specific effect estimates for the primary outcome",
  9.0
)
ft_or11_sens <- make_ft(
  or11_sensitivity_results %>%
    transmute(Analysis = label,
              `Clones, n` = fmt_n(n_clone),
              `Events, n` = fmt_n(n_event),
              `HR (95% CI)` = HR_95CI),
  "Online Resource 11. Sensitivity of the primary outcome to database and alternative analytic definitions",
  8.3
)
ft_or12_postvp <- ft_post_vp_events_bundle
ft_or12_lactate <- ft_lactate_traj
ft_or13_nee <- make_ft(
  nee_traj_table,
  "Online Resource 13. Weighted norepinephrine-equivalent dose trajectory from time zero through 72 hours",
  8.2
)

save_or_figure(fig_or7_vp_timing, "OR7_vp_timing_distribution")
save_or_figure(fig_or8_weights, "OR8_ipcw_weight_distribution")
save_or_figure(fig_or10_database, "OR10_database_specific_forest")
save_or_figure(fig_or11_sensitivity, "OR11_sensitivity_forest")
save_or_figure(fig_lactate_traj, "OR12_lactate_trajectory_after_vp")
save_or_figure(fig_or13_nee_traj, "OR13_weighted_nee_trajectory")
save_or_figure(fig_or15_baseline_dist, "OR15_baseline_severity_distributions")
save_or_figure(fig_or16_hemo_traj, "OR16_hemodynamic_trajectories", height = 8.2)

write_csv(vp_timing_bins, file.path(OR_REVISED_DIR, "OR7_vp_timing_distribution.csv"))
write_csv(or8_weight_label, file.path(OR_REVISED_DIR, "OR8_ipcw_weight_summary.csv"))
write_csv(database_effects, file.path(OR_REVISED_DIR, "OR10_database_specific_estimates.csv"))
write_csv(or11_sensitivity_results, file.path(OR_REVISED_DIR, "OR11_sensitivity_estimates.csv"))
write_csv(lactate_trajectory_table, file.path(OR_REVISED_DIR, "OR12_lactate_trajectory_after_vp.csv"))
write_csv(post_vp_events_bundle_table, file.path(OR_REVISED_DIR, "OR12_post_vp_events.csv"))
write_csv(nee_traj_table, file.path(OR_REVISED_DIR, "OR13_weighted_nee_trajectory_table.csv"))
write_csv(or14_table, file.path(OR_REVISED_DIR, "OR14_catecholamine_vfd_table.csv"))
write_csv(or15_dist_df, file.path(OR_REVISED_DIR, "OR15_baseline_distribution_plot_data.csv"))
write_csv(hemo_traj_summary, file.path(OR_REVISED_DIR, "OR16_hemodynamic_trajectory_summary.csv"))

or7_16_docx <- file.path(OR_REVISED_DIR, "VP_online_resources_7_16_revised_only.docx")
save_as_docx(
  `Online Resource 10 database estimates` = ft_or10_db,
  `Online Resource 11 sensitivity estimates` = ft_or11_sens,
  `Online Resource 12 lactate trajectory` = ft_or12_lactate,
  `Online Resource 12 post-VP events` = ft_or12_postvp,
  `Online Resource 13 NEE trajectory table` = ft_or13_nee,
  `Online Resource 14 catecholamine and VFD` = ft_or14,
  path = or7_16_docx,
  pr_section = prop_section(
    page_size = page_size(orient = "landscape"),
    page_margins = page_mar(top = 0.45, bottom = 0.45, left = 0.45, right = 0.45)
  )
)

or7_16_pptx <- file.path(OR_REVISED_DIR, "VP_online_resources_7_16_figures_editable.pptx")
or7_16_ppt <- read_pptx()
or7_16_ppt <- add_figure_slide(or7_16_ppt, "Online Resource 7. Time to first vasopressin initiation", fig_or7_vp_timing)
or7_16_ppt <- add_figure_slide(or7_16_ppt, "Online Resource 8. IPCW distribution", fig_or8_weights)
or7_16_ppt <- add_figure_slide(or7_16_ppt, "Online Resource 10. Database-specific estimates", fig_or10_database)
or7_16_ppt <- add_figure_slide(or7_16_ppt, "Online Resource 11. Sensitivity analyses", fig_or11_sensitivity)
or7_16_ppt <- add_figure_slide(or7_16_ppt, "Online Resource 12. Lactate trajectory after VP", fig_lactate_traj)
or7_16_ppt <- add_figure_slide(or7_16_ppt, "Online Resource 13. Weighted NEE trajectory", fig_or13_nee_traj)
or7_16_ppt <- add_figure_slide(or7_16_ppt, "Online Resource 15. Baseline severity distributions", fig_or15_baseline_dist)
or7_16_ppt <- add_figure_slide(or7_16_ppt, "Online Resource 16. Hemodynamic trajectories", fig_or16_hemo_traj, h = 7.2)
print(or7_16_ppt, target = or7_16_pptx)

or7_16_index <- tibble(
  file = c(
    "OR7_vp_timing_distribution.png",
    "OR8_ipcw_weight_distribution.png",
    "OR10_database_specific_forest.png",
    "OR11_sensitivity_forest.png",
    "OR12_lactate_trajectory_after_vp.png",
    "OR13_weighted_nee_trajectory.png",
    "OR15_baseline_severity_distributions.png",
    "OR16_hemodynamic_trajectories.png",
    basename(or7_16_docx),
    basename(or7_16_pptx)
  ),
  description = c(
    "Revised Online Resource 7 figure",
    "Revised Online Resource 8 figure",
    "Revised Online Resource 10 figure",
    "Revised Online Resource 11 figure with database-specific and time-zero sensitivity estimates",
    "Revised Online Resource 12 post-VP lactate figure",
    "Revised Online Resource 13 figure",
    "Revised Online Resource 15 figure",
    "Revised Online Resource 16 figure",
    "Editable Word file containing revised OR7-16 tables",
    "Editable PowerPoint file containing revised OR7-16 figures"
  )
)
write_csv(or7_16_index, file.path(OR_REVISED_DIR, "OR7_16_output_index.csv"))

make_definition_ft <- function(dat, title, widths, font_size = 6.6) {
  ft <- flextable(dat)
  ft <- theme_vanilla(ft)
  ft <- font(ft, fontname = FONT_FAMILY, part = "all")
  ft <- fontsize(ft, size = font_size, part = "all")
  ft <- fontsize(ft, size = font_size + 0.8, part = "header")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "left", part = "all")
  ft <- valign(ft, valign = "top", part = "all")
  ft <- padding(ft, padding.top = 3, padding.bottom = 3,
                padding.left = 4, padding.right = 4, part = "all")
  for (j in seq_along(widths)) ft <- width(ft, j = j, width = widths[[j]])
  ft <- add_header_lines(ft, values = title)
  ft <- bold(ft, i = 1, part = "header")
  ft <- fontsize(ft, i = 1, size = font_size + 2, part = "header")
  ft <- set_table_properties(ft, layout = "fixed", width = 1)
  ft
}

definition_comorbidity_safety <- tribble(
  ~Domain, ~Concept, ~Source, ~Definition,
  "Comorbidity", "Chronic cardiac disease", "MIMIC-IV",
  "ICD-9: 39891, 402*, 404*, 4254-4259, 428*, 412*, 414*. ICD-10: I50*, I110, I130, I132, I42*, I43*, I25*.",
  "Comorbidity", "Chronic cardiac disease", "eICU-CRD",
  "Medical-history text terms: heart failure, cardiomyopathy, coronary artery, ischemic/ischaemic heart, CHF, HFrEF, HFpEF, or valvular disease.",
  "Comorbidity", "Hypertension", "MIMIC-IV",
  "ICD-9: 401*, 402*, 403*, 404*, 405*. ICD-10: I10*, I11*, I12*, I13*, I15*.",
  "Comorbidity", "Hypertension", "eICU-CRD",
  "Medical-history text terms: hypertension or HTN.",
  "Comorbidity", "Chronic kidney disease", "MIMIC-IV",
  "ICD-9: 585*, V420, V451, V56*. ICD-10: N18*, Z940, Z992, Z49*.",
  "Comorbidity", "Chronic kidney disease", "eICU-CRD",
  "Medical-history text terms: chronic kidney disease, CKD, end-stage renal disease, ESRD, or dialysis dependent.",
  "Comorbidity", "Diabetes mellitus", "MIMIC-IV",
  "ICD-9: 250*. ICD-10: E10*, E11*, E12*, E13*, E14*.",
  "Comorbidity", "Diabetes mellitus", "eICU-CRD",
  "Medical-history text terms: diabetes, DM type, NIDDM, IDDM, type 1 DM, type 2 DM, T1DM, or T2DM.",
  "Comorbidity", "Malignancy", "MIMIC-IV",
  "ICD-9: 140*-165* and 170*-208*. ICD-10: C00*-C99*.",
  "Comorbidity", "Malignancy", "eICU-CRD",
  "Medical-history text terms: cancer, malignancy, carcinoma, lymphoma, leukemia, metastasis/metastatic, sarcoma, or myeloma.",
  "Comorbidity", "HIV/AIDS", "MIMIC-IV",
  "ICD-9: 042, 043, 044, V08. ICD-10: B20*, B21*, B22*, B23*, B24*, Z21.",
  "Comorbidity", "HIV/AIDS", "eICU-CRD",
  "Medical-history text terms: HIV or AIDS.",
  "Comorbidity", "Dementia", "MIMIC-IV",
  "ICD-9: 290*, 2941, 3310, 3311, 3312, 33182. ICD-10: F00*, F01*, F02*, F03*, G30*, G310, G311, G3183.",
  "Comorbidity", "Dementia", "eICU-CRD",
  "Medical-history text terms: dementia, Alzheimer disease, or cognitive impairment.",
  "Exploratory safety", "Deep vein thrombosis", "MIMIC-IV",
  "ICD-9: 4511*, 4512*, 4532*, 4534*, 4538*, 4539*. ICD-10: I801*, I802*, I824*, I825*, I826*, I828*, I829*.",
  "Exploratory safety", "Deep vein thrombosis", "eICU-CRD",
  "Diagnosis text terms: deep vein thrombosis, DVT, or venous thrombosis; ICD codes same as MIMIC-IV when available.",
  "Exploratory safety", "Pulmonary embolism", "MIMIC-IV",
  "ICD-9: 4151*. ICD-10: I26*.",
  "Exploratory safety", "Pulmonary embolism", "eICU-CRD",
  "Diagnosis text terms: pulmonary embolism or pulmonary emboli; ICD codes same as MIMIC-IV when available.",
  "Exploratory safety", "Mesenteric ischemia", "MIMIC-IV",
  "ICD-9: 5570, 5571, 5579. ICD-10: K550*, K551*, K559*.",
  "Exploratory safety", "Mesenteric ischemia", "eICU-CRD",
  "Diagnosis text terms: mesenteric ischemia, bowel ischemia, intestinal ischemia, or acute mesenteric ischemia; ICD codes same as MIMIC-IV when available.",
  "Exploratory safety", "Digital/peripheral ischemia", "MIMIC-IV",
  "ICD-9: 7854, 44024, 4439. ICD-10: I96, I9981, I7024*, I7025*, I7026*, I739.",
  "Exploratory safety", "Digital/peripheral ischemia", "eICU-CRD",
  "Diagnosis text terms: digital ischemia, peripheral ischemia, limb ischemia, finger ischemia, toe ischemia, gangrene, or necrosis; ICD codes same as MIMIC-IV when available."
)

definition_outcomes_treatments <- tribble(
  ~Domain, ~Concept, ~Definition, ~Notes,
  "Eligibility/time zero", "Primary time zero",
  "First hour at which the norepinephrine infusion rate reached at least 0.25 mcg/kg/min, with at least 3 qualifying hours within the subsequent 3-hour window and no prior vasopressin.",
  "Used for the revised main analysis.",
  "Eligibility/time zero", "NEE time-zero sensitivity",
  "First hour at which the norepinephrine-equivalent dose reached at least 0.25 mcg/kg/min, using the same sustained-window rule.",
  "Used as a sensitivity analysis after moving the main analysis to norepinephrine dose.",
  "Eligibility/time zero", "Norepinephrine-only time-zero sensitivity",
  "First hour at which norepinephrine reached at least 0.25 mcg/kg/min while epinephrine, dopamine, and phenylephrine were not running, again requiring at least 3 qualifying hours within the subsequent 3-hour window.",
  "Used to address the concern that non-norepinephrine first-line vasopressors may deviate from guideline-concordant practice.",
  "Treatment strategy", "Ultra-early VP strategy",
  "Vasopressin initiation within 0-3 h after time zero.",
  "Patients were cloned into strategy arms and artificially censored at nonadherence.",
  "Treatment strategy", "Early VP strategy",
  "Vasopressin initiation within >3-6 h after time zero.",
  "Patients receiving vasopressin after 6 h or not receiving vasopressin during the grace period were censored for the early strategy at 6 h.",
  "Primary outcome", "28-day mortality",
  "Death within 28 days after time zero.",
  "Primary effect estimated with clone-censor-weight Cox models; risk difference and RMST difference were also derived.",
  "Secondary outcome", "AKI stage >=2 or >=3",
  "First time after time zero at which KDIGO stage >=2 or >=3 was met based on creatinine and urine-output criteria.",
  "Patients already meeting the relevant AKI stage at time zero were excluded from that secondary-outcome risk set.",
  "Secondary outcome", "Renal replacement therapy",
  "First initiation after time zero of RRT/CRRT. MIMIC-IV used procedure/input events including dialysis/CRRT-related item IDs; eICU-CRD used treatment/intake-output terms including RRT, dialysis, ultrafiltration, CAVHD, CVVH, or SLED.",
  "Patients with pre-time-zero RRT/CRRT were excluded from that secondary-outcome risk set.",
  "Secondary outcome", "Medically treated arrhythmia",
  "First initiation after time zero of an antiarrhythmic treatment, operationalized from medication records for amiodarone, lidocaine, or procainamide in the source derivation.",
  "Patients with pre-time-zero antiarrhythmic treatment were excluded from that secondary-outcome risk set.",
  "Secondary outcome", "Shock resolution",
  "First time after time zero at which vasopressor-free status was maintained for at least 24 h.",
  "Estimated as a time-to-event secondary outcome.",
  "Secondary outcome", "Net negative fluid balance",
  "First time after time zero at which cumulative net fluid balance became negative.",
  "Estimated as a time-to-event secondary outcome.",
  "Baseline/treatment", "Fluid input and balance",
  "Hourly fluid input and output were summarized in the 24 h and 6 h before time zero, the 24 h and 6 h before VP initiation, and the 6 h after VP initiation.",
  "MIMIC-IV used inputevents/outputevents; eICU-CRD used intake-output records.",
  "Baseline/treatment", "Steroid timing",
  "Corticosteroid exposure was summarized before time zero, before VP initiation, and after VP initiation using available hourly treatment indicators or original source flags when hourly corticosteroid data were unavailable.",
  "Used for baseline and descriptive treatment summaries.",
  "VP exposure", "Initial VP dose",
  "First recorded vasopressin infusion rate at VP initiation.",
  "Also categorized as <0.03, approximately 0.03, or >0.03 U/min.",
  "VP exposure", "Maximum and cumulative VP dose",
  "Maximum and cumulative vasopressin dose during the first 6 h and first 24 h after VP initiation.",
  "Reported by observed timing group.",
  "Post-VP signal", "Early VP interruption",
  "Any recorded vasopressin infusion rate of 0 during 0-3 h after VP initiation.",
  "Compared between adherent ultra-early and early strategy clones using IPCW-weighted percentages.",
  "Post-VP signal", "Permanent VP discontinuation within 3 h",
  "First vasopressin off-hour within 3 h after VP initiation followed by no restart through 24 h.",
  "Compared between adherent ultra-early and early strategy clones using IPCW-weighted percentages.",
  "Post-VP signal", "New dopamine/dobutamine/epinephrine",
  "No use during the 24 h before or at VP initiation, followed by new use within 3 h, 6 h, or 24 h after VP initiation.",
  "Compared between adherent ultra-early and early strategy clones using IPCW-weighted percentages.",
  "Post-VP signal", "Lactate trajectory after VP",
  "Lactate values at VP initiation and 3, 6, 12, and 24 h after VP initiation.",
  "Displayed as IPCW-weighted median and IQR among adherent ultra-early and early strategy clones.",
  "Exploratory safety", "Safety model estimates",
  "Diagnosis-based DVT, PE, mesenteric ischemia, and digital/peripheral ischemia were compared using IPCW-weighted Poisson models with robust clustered standard errors.",
  "HRs were not estimated because the current safety definitions are diagnosis based and do not contain reliable event onset times."
)

ft_definition_comorbidity_safety <- make_definition_ft(
  definition_comorbidity_safety,
  "Supplementary Table. Definitions of comorbidities and exploratory safety outcomes",
  widths = c(1.05, 1.75, 0.95, 6.25),
  font_size = 6.4
)
ft_definition_outcomes_treatments <- make_definition_ft(
  definition_outcomes_treatments,
  "Supplementary Table. Definitions of time zero, outcomes, and treatment-related variables",
  widths = c(1.15, 1.95, 4.65, 2.25),
  font_size = 6.4
)

supp_definitions_docx <- file.path(SUPP_BUNDLE_DIR, "supplementary_definitions_editable.docx")
save_as_docx(
  `Comorbidities and safety definitions` = ft_definition_comorbidity_safety,
  `Outcome and treatment definitions` = ft_definition_outcomes_treatments,
  path = supp_definitions_docx,
  pr_section = prop_section(
    page_size = page_size(orient = "landscape"),
    page_margins = page_mar(top = 0.45, bottom = 0.45, left = 0.45, right = 0.45)
  )
)

supp_caption_entries <- tribble(
  ~Item, ~Caption, ~Abbreviations,
  "Supplementary Table: Baseline ultra-early, early, and other",
  "Baseline characteristics according to observed VP timing: ultra-early, early, and other. Other includes VP initiation after 6 h and no observed VP initiation.",
  "IQR, interquartile range; NE, norepinephrine; VP, vasopressin.",
  "Supplementary Table/Figure: Sensitivity analyses",
  "Sensitivity analyses for 28-day mortality after changing the time-zero definition to NEE, restricting to Nad-only at time zero, and excluding patients receiving dobutamine at time zero.",
  "CI, confidence interval; HR, hazard ratio; IPCW, inverse-probability-of-censoring weighting; Nad, norepinephrine; NEE, norepinephrine-equivalent dose.",
  "Supplementary Table: Vasopressin dose",
  "Initial, maximum, and cumulative vasopressin dose according to observed VP timing.",
  "IQR, interquartile range; VP, vasopressin.",
  "Supplementary Table: Survivors versus non-survivors",
  "Baseline characteristics of 28-day survivors and non-survivors.",
  "IQR, interquartile range; MAP, mean arterial pressure; NE, norepinephrine; SOFA, Sequential Organ Failure Assessment.",
  "Supplementary Table/Figure: Post-VP physiologic and cardiovascular support signals",
  "IPCW-weighted comparison of lactate trajectory, early VP interruption, and new dopamine, dobutamine, or epinephrine initiation after VP initiation among adherent ultra-early and early strategy clones.",
  "IPCW, inverse-probability-of-censoring weighting; IQR, interquartile range; VP, vasopressin.",
  "Supplementary Table/Figure: Safety outcomes",
  "Exploratory diagnosis-based safety outcomes, including DVT, PE, mesenteric ischemia, and digital/peripheral ischemia. Descriptive results are shown by observed timing; model estimates compare the ultra-early and early strategies using risk ratios rather than hazard ratios because event onset times are unavailable.",
  "CI, confidence interval; DVT, deep vein thrombosis; PE, pulmonary embolism; RR, risk ratio; VP, vasopressin.",
  "Supplementary Table: Definitions",
  "Definitions of comorbidities, secondary outcomes, treatment-related variables, and exploratory safety outcomes, including ICD codes and eICU text patterns where applicable.",
  "AKI, acute kidney injury; CRRT, continuous renal replacement therapy; ICD, International Classification of Diseases; RRT, renal replacement therapy."
)
ft_supp_captions <- make_caption_ft(supp_caption_entries, 7.7)

supp_tables_docx <- file.path(SUPP_BUNDLE_DIR, "supplementary_tables_editable.docx")
save_as_docx(
  `Baseline ultra-early early other` = ft_baseline_three,
  `Sensitivity analyses` = ft_requested_sens,
  `Vasopressin dose` = ft_vp_dose_bundle,
  `Survivor versus non-survivor` = ft_survivor_bundle,
  `Lactate trajectory` = ft_lactate_traj,
  `Post-VP events` = ft_post_vp_events_bundle,
  `Safety descriptive` = ft_safety_desc_bundle,
  `Safety model` = ft_safety_model,
  `Comorbidities and safety definitions` = ft_definition_comorbidity_safety,
  `Outcome and treatment definitions` = ft_definition_outcomes_treatments,
  `Captions and abbreviations` = ft_supp_captions,
  path = supp_tables_docx,
  pr_section = prop_section(
    page_size = page_size(orient = "landscape"),
    page_margins = page_mar(top = 0.45, bottom = 0.45, left = 0.45, right = 0.45)
  )
)

supp_figures_pptx <- file.path(SUPP_BUNDLE_DIR, "supplementary_figures_editable.pptx")
supp_bundle_ppt <- read_pptx()
supp_bundle_ppt <- add_figure_slide(supp_bundle_ppt, "Supplementary Figure. Requested sensitivity analyses", fig_requested_sens)
supp_bundle_ppt <- add_figure_slide(supp_bundle_ppt, "Supplementary Figure. Lactate trajectory after vasopressin", fig_lactate_traj)
supp_bundle_ppt <- add_figure_slide(supp_bundle_ppt, "Supplementary Figure. Exploratory safety outcome models", fig_safety_model)
print(supp_bundle_ppt, target = supp_figures_pptx)

ggsave(file.path(SUPP_BUNDLE_DIR, "figure_sensitivity_timezero_dobutamine.pdf"),
       fig_requested_sens, width = 11, height = 7.5, device = cairo_pdf, bg = "white")
ggsave(file.path(SUPP_BUNDLE_DIR, "figure_lactate_trajectory_after_vp.pdf"),
       fig_lactate_traj, width = 11, height = 7.5, device = cairo_pdf, bg = "white")
ggsave(file.path(SUPP_BUNDLE_DIR, "figure_safety_outcome_model.pdf"),
       fig_safety_model, width = 11, height = 7.5, device = cairo_pdf, bg = "white")
ggsave(file.path(SUPP_BUNDLE_DIR, "figure_sensitivity_timezero_dobutamine.png"),
       fig_requested_sens, width = 11, height = 7.5, dpi = 220, bg = "white")
ggsave(file.path(SUPP_BUNDLE_DIR, "figure_lactate_trajectory_after_vp.png"),
       fig_lactate_traj, width = 11, height = 7.5, dpi = 220, bg = "white")
ggsave(file.path(SUPP_BUNDLE_DIR, "figure_safety_outcome_model.png"),
       fig_safety_model, width = 11, height = 7.5, dpi = 220, bg = "white")

write_csv(table1_three_body %>% select(-.section), file.path(SUPP_BUNDLE_DIR, "table_baseline_ultra_early_other.csv"))
write_csv(requested_sens_table, file.path(SUPP_BUNDLE_DIR, "table_sensitivity_timezero_dobutamine.csv"))
write_csv(vp_dose_bundle_table, file.path(SUPP_BUNDLE_DIR, "table_vasopressin_dose.csv"))
write_csv(survivor_bundle_table, file.path(SUPP_BUNDLE_DIR, "table_survivor_non_survivor.csv"))
write_csv(lactate_trajectory_table, file.path(SUPP_BUNDLE_DIR, "table_lactate_trajectory_after_vp.csv"))
write_csv(post_vp_events_bundle_table, file.path(SUPP_BUNDLE_DIR, "table_post_vp_interruption_new_agents.csv"))
write_csv(safety_desc_bundle_table, file.path(SUPP_BUNDLE_DIR, "table_safety_descriptive_by_timing.csv"))
write_csv(safety_model_table, file.path(SUPP_BUNDLE_DIR, "table_safety_model_estimates.csv"))
write_csv(definition_comorbidity_safety, file.path(SUPP_BUNDLE_DIR, "table_definitions_comorbidities_safety.csv"))
write_csv(definition_outcomes_treatments, file.path(SUPP_BUNDLE_DIR, "table_definitions_outcomes_treatments.csv"))
write_csv(supp_caption_entries, file.path(SUPP_BUNDLE_DIR, "captions_abbreviations.csv"))

supp_output_index <- tibble(
  file = c(
    basename(supp_tables_docx),
    basename(supp_definitions_docx),
    basename(supp_figures_pptx),
    "table_baseline_ultra_early_other.csv",
    "table_sensitivity_timezero_dobutamine.csv",
    "table_vasopressin_dose.csv",
    "table_survivor_non_survivor.csv",
    "table_lactate_trajectory_after_vp.csv",
    "table_post_vp_interruption_new_agents.csv",
    "table_safety_descriptive_by_timing.csv",
    "table_safety_model_estimates.csv",
    "table_definitions_comorbidities_safety.csv",
    "table_definitions_outcomes_treatments.csv",
    "figure_sensitivity_timezero_dobutamine.png",
    "figure_lactate_trajectory_after_vp.png",
    "figure_safety_outcome_model.png",
    "captions_abbreviations.csv"
  ),
  description = c(
    "Editable Word file containing all requested supplementary tables and captions",
    "Editable Word file containing definitions of comorbidities, outcomes, treatments, and exploratory safety outcomes",
    "Editable PowerPoint file containing all requested supplementary figures",
    "Baseline comparison: ultra-early, early, and other observed groups",
    "Sensitivity analyses: NEE time zero, Nad-only time zero, and dobutamine exclusion",
    "Vasopressin dose summaries by observed timing",
    "Baseline comparison of 28-day survivors and non-survivors",
    "IPCW-weighted lactate trajectory after vasopressin initiation among adherent assigned-strategy clones",
    "IPCW-weighted early VP interruption and new dopamine/dobutamine/epinephrine after VP",
    "Exploratory safety outcomes by observed timing",
    "Exploratory weighted model estimates for safety outcomes",
    "Definitions of comorbidities and exploratory safety outcomes, including ICD codes and eICU text patterns",
    "Definitions of time zero, secondary outcomes, treatment exposures, and post-VP signals",
    "PNG preview of sensitivity forest plot",
    "PNG preview of lactate trajectory plot",
    "PNG preview of safety model forest plot",
    "Captions and abbreviations for the requested supplementary items"
  )
)
write_csv(supp_output_index, file.path(SUPP_BUNDLE_DIR, "output_index.csv"))
supp_bundle_expected <- file.path(SUPP_BUNDLE_DIR, supp_output_index$file)

qa_summary <- tibble(
  item = c(
    "N eligible",
    "MIMIC-IV eligible",
    "eICU-CRD eligible",
    "Adherent ultra 0-3 h",
    "Adherent early >3-6 h",
    "Never VP within 6 h",
    "Primary HR",
    "Primary weighted risk ultra",
    "Primary weighted risk early",
    "Primary RD pp",
    "Primary RMST diff days"
  ),
  value = c(
    fmt_n(nrow(pt_ne)),
    fmt_n(sum(pt_ne$src == "miiv")),
    fmt_n(sum(pt_ne$src == "eicu")),
    fmt_n(sum(!is.na(pt_ne$tvp_rel) & pt_ne$tvp_rel <= 3)),
    fmt_n(sum(!is.na(pt_ne$tvp_rel) & pt_ne$tvp_rel > 3 & pt_ne$tvp_rel <= 6)),
    fmt_n(sum(is.na(pt_ne$tvp_rel) | pt_ne$tvp_rel > 6)),
    primary_sens$HR_95CI[primary_sens$label == "Primary: NE time zero"],
    sprintf("%.3f", pt_rd$r_u),
    sprintf("%.3f", pt_rd$r_e),
    sprintf("%.2f (%.2f to %.2f)", pt_rd$rd * 100, ci_rd[1] * 100, ci_rd[2] * 100),
    sprintf("%.2f (%.2f to %.2f)", pt_rmst, ci_rmst[1], ci_rmst[2])
  )
)
write_csv(qa_summary, file.path(OUT_DIR, "summary_ne_main_outputs.csv"))

expected <- c(
  fig_pptx, supp_fig_pptx, tables_docx, online_docx, captions_docx,
  file.path(TAB_DIR, "table1_ne_main.docx"),
  file.path(TAB_DIR, "table_primary_sensitivity_ne_main.docx"),
  file.path(TAB_DIR, "table_secondary_ne_main.docx"),
  file.path(TAB_DIR, "supplementary_tables_ne_main.docx"),
  file.path(FIG_DIR, "figure1_flow_ne_main.pdf"),
  file.path(FIG_DIR, "figure2_cuminc_ne_main.pdf"),
  file.path(FIG_DIR, "figure3_outcomes_ne_main.pdf"),
  file.path(FIG_DIR, "figure4_subgroup_ne_main.pdf"),
  file.path(FIG_DIR, "supp_figure_primary_sensitivity_ne_main.pdf"),
  file.path(FIG_DIR, "supp_figure1_vp_timing_ne_main.pdf"),
  file.path(FIG_DIR, "supp_figure2_ipcw_weights_ne_main.pdf"),
  file.path(FIG_DIR, "supp_figure3_uncensored_clones_ne_main.pdf"),
  file.path(QA_DIR, "figure1_flow_ne_main.png"),
  file.path(QA_DIR, "figure2_cuminc_ne_main.png"),
  file.path(QA_DIR, "figure3_outcomes_ne_main.png"),
  file.path(QA_DIR, "figure4_subgroup_ne_main.png"),
  file.path(QA_DIR, "supp_figure_primary_sensitivity_ne_main.png"),
  file.path(QA_DIR, "supp_figure1_vp_timing_ne_main.png"),
  file.path(QA_DIR, "supp_figure2_ipcw_weights_ne_main.png"),
  file.path(QA_DIR, "supp_figure3_uncensored_clones_ne_main.png"),
  supp_bundle_expected
)
missing <- expected[!file.exists(expected) | file.info(expected)$size <= 0]
if (length(missing)) {
  stop("Missing or empty outputs:\n", paste(missing, collapse = "\n"))
}
if (any(str_detect(primary_sens$label, regex("no time-varying", ignore_case = TRUE)))) {
  stop("An excluded IPCW sensitivity unexpectedly appeared in primary/sensitivity outputs.")
}

message("Done. Outputs written to: ", OUT_DIR)
print(qa_summary)
