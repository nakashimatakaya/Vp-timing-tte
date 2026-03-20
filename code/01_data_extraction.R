################################################################################
# 01_data_extraction.R
#
# Cohort extraction and hourly harmonization using ricu.
# Sources: MIMIC-IV (miiv), eICU-CRD (eicu)
#
# Outputs:
#   - septic_shock_imputed_df_new (hourly imputed time-series)
#   - vp_ts_out_secondary_*.rds  (uid-level outcome/subgroup table)
#   - vp_pt_vp6_*.rds            (patient-level analytic cohort)
################################################################################

library(ricu)
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(data.table)
library(purrr)
library(slider)
library(stringr)

# ==============================================================================
# 0. Paths (update to your local environment)
# ==============================================================================

DATA_DIR  <- "/Users/nakashimatakaya/Desktop/SID/Data/Analysis_Data"
CACHE_DIR <- file.path(DATA_DIR, "_cache_vp_secondary")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

TS_OUT_PATH <- file.path(DATA_DIR, "vp_ts_out_secondary_20260223.rds")

# ==============================================================================
# 1. Attach data sources
# ==============================================================================

src_data_avail()
attach_src("miiv")
attach_src("eicu")

# ==============================================================================
# 2. Variable definitions
# ==============================================================================

demo_vars    <- c("age", "sex", "weight", "height", "bmi")
vaso_vars    <- c("norepi_rate", "norepi_equiv", "adh_rate", "vaso_ind")
vaso_addon   <- c("epi_rate", "dopa_rate", "dobu_rate", "phn_rate")
vital_vars   <- c("map", "sbp", "dbp", "hr")
lab_vars     <- c("na", "k", "cl", "mg", "phos", "ca", "cai",
                   "lact", "alb", "crea", "ph", "bicar", "tco2",
                   "pco2", "be")
score_vars   <- c("sofa", "sofa_cardio", "sep3", "abx", "cort",
                   "vent_ind", "mech_vent", "urine24")
outcome_vars <- c("death", "los_icu", "los_hosp")

vars_all <- unique(c(demo_vars, vaso_vars, vital_vars, lab_vars,
                     score_vars, outcome_vars))

# ==============================================================================
# 3. Helper functions for ricu output standardization
# ==============================================================================

get_id_col <- function(x) {
  v <- tryCatch(ricu::id_vars(x), error = function(e) character(0))
  if (length(v) == 0) return(NA_character_)
  v[length(v)]
}

get_time_col <- function(x) {
  is_ts <- tryCatch(ricu::is_ts_tbl(x), error = function(e) FALSE)
  if (!is_ts) return(NA_character_)
  tryCatch(ricu::index_var(x), error = function(e) NA_character_)
}

to_hours_safe <- function(t, src) {
  if (inherits(t, "difftime")) return(as.numeric(t, units = "hours"))
  if (is.numeric(t)) {
    if (src == "eicu" && suppressWarnings(max(t, na.rm = TRUE)) > 500) return(t / 60)
    return(t)
  }
  suppressWarnings(as.numeric(t))
}

std_tbl <- function(x, src) {
  if (is.null(x)) return(tibble())
  id_col <- get_id_col(x)
  t_col  <- get_time_col(x)
  df <- as_tibble(x)
  if (!is.na(id_col) && id_col %in% names(df)) {
    df <- df %>% rename(all_id = all_of(id_col))
  } else {
    df$all_id <- NA
  }
  if (!is.na(t_col) && t_col %in% names(df)) {
    df <- df %>%
      rename(time_raw = all_of(t_col)) %>%
      mutate(time_hr = to_hours_safe(time_raw, src)) %>%
      select(-time_raw)
  } else {
    df$time_hr <- NA_real_
  }
  df %>% mutate(src = src)
}

ensure_cols <- function(df, cols) {
  for (c in cols) {
    if (!c %in% names(df)) df[[c]] <- NA
  }
  df
}

# ==============================================================================
# 4. Load data from each source (uncomment to run from scratch)
# ==============================================================================

# --- MIMIC-IV ---
# miiv_raw <- load_concepts(vars_all, "miiv", merge_data = TRUE, verbose = FALSE)
# miiv_raw <- tryCatch(fill_gaps(miiv_raw), error = function(e) miiv_raw)
# miiv_df  <- std_tbl(miiv_raw, "miiv") %>% ensure_cols(vars_all)
# write_rds(miiv_df, file.path(DATA_DIR, "miiv_df"))

# --- eICU-CRD ---
# eicu_raw <- load_concepts(vars_all, "eicu", merge_data = TRUE, verbose = FALSE)
# eicu_raw <- tryCatch(fill_gaps(eicu_raw), error = function(e) eicu_raw)
# eicu_df  <- std_tbl(eicu_raw, "eicu") %>% ensure_cols(vars_all)
# write_rds(eicu_df, file.path(DATA_DIR, "eicu_df"))

# --- Additional vasopressor data ---
# vaso_addiction <- load_concepts(vaso_addon, c("miiv", "eicu"),
#                                 merge_data = TRUE, verbose = FALSE)
# write_rds(vaso_addiction, file.path(DATA_DIR, "vaso_addiction"))

# ==============================================================================
# 5. Load cached data and standardize
# ==============================================================================

miiv_df  <- read_rds(file.path(DATA_DIR, "miiv_df"))
eicu_df3 <- read_rds(file.path(DATA_DIR, "eicu_df3"))
vaso_addiction <- read_rds(file.path(DATA_DIR, "vaso_addiction"))

align_df <- function(df) {
  df %>%
    mutate(all_id = as.character(all_id)) %>%
    ensure_cols(vars_all) %>%
    select(all_id, all_of(vars_all), time_hr, src)
}

miiv_df <- align_df(miiv_df)
eicu_df <- align_df(eicu_df3)

# ==============================================================================
# 6. Filter to septic shock cohort (sep3 + NE)
# ==============================================================================

filter_ts_sep3_ne <- function(df) {
  df2 <- df %>%
    mutate(
      sep3  = as.logical(sep3),
      ne_on = (!is.na(norepi_rate) & norepi_rate > 0) |
              (!is.na(norepi_equiv) & norepi_equiv > 0)
    )
  ids_keep <- df2 %>%
    group_by(src, all_id) %>%
    summarise(has_sep3 = any(sep3, na.rm = TRUE),
              has_ne = any(ne_on, na.rm = TRUE),
              .groups = "drop") %>%
    filter(has_sep3 & has_ne) %>%
    select(src, all_id)
  df2 %>%
    semi_join(ids_keep, by = c("src", "all_id")) %>%
    select(-ne_on)
}

miiv_ts_keep <- filter_ts_sep3_ne(miiv_df)
eicu_ts_keep <- filter_ts_sep3_ne(eicu_df)

# ==============================================================================
# 7. Pool databases
# ==============================================================================

ts_all <- bind_rows(miiv_ts_keep, eicu_ts_keep) %>%
  filter(!is.na(all_id)) %>%
  left_join(
    vaso_addiction %>%
      mutate(stay_id = as.character(stay_id),
             starttime = as.numeric(starttime)),
    by = c("all_id" = "stay_id", "src" = "source", "time_hr" = "starttime")
  )

ts_all_ts <- ts_all %>%
  filter(!is.na(time_hr)) %>%
  arrange(src, all_id, time_hr)

ts_all_ts %>%
  group_by(src) %>%
  summarise(n_rows = n(), n_ids = n_distinct(paste(src, all_id)))

# ==============================================================================
# 8. Imputation parameters
# ==============================================================================

PARAM <- list(
  grace_hr        = 6,
  max_followup_hr = 28 * 24,
  scd_window_hr   = 2,
  cap_vitals_hr   = 2,
  cap_fastlab_hr  = 6,
  cap_slowlab_hr  = 24,
  cap_sofa_hr     = 24,
  eps_rate         = 0
)

# ==============================================================================
# 9. LOCF imputation with cap
# ==============================================================================

id_cols  <- c("src", "all_id")
time_col <- "time_hr"

rate_vars    <- intersect(c("norepi_rate", "norepi_equiv", "adh_rate",
                            "epi_rate", "dopa_rate", "phn_rate", "dobu_rate"),
                          names(ts_all_ts))
vitals_vars  <- intersect(c("map", "sbp", "dbp", "hr"), names(ts_all_ts))
fastlab_vars <- intersect(c("lact", "ph", "pco2", "be", "bicar", "tco2"), names(ts_all_ts))
slowlab_vars <- intersect(c("na", "k", "cl", "mg", "phos", "ca", "cai",
                            "crea", "alb", "urine24"), names(ts_all_ts))
sofa_vars    <- intersect(c("sofa", "sofa_cardio"), names(ts_all_ts))
static_num   <- intersect(c("age", "weight", "height", "bmi", "los_icu", "los_hosp"),
                           names(ts_all_ts))
flag_vars    <- intersect(c("sep3", "abx", "cort", "vent_ind", "vaso_ind", "death"),
                          names(ts_all_ts))

impute_locf_cap_keepna <- function(dt, var, cap_hr,
                                   id_cols = c("src", "all_id"),
                                   time_col = "time_hr") {
  if (!var %in% names(dt)) return(invisible(NULL))
  suppressWarnings(dt[, (var) := as.numeric(get(var))])
  tlast <- paste0(var, "__tlast")
  dt[, (tlast) := fifelse(!is.na(get(var)), get(time_col), as.numeric(NA)), by = id_cols]
  dt[, (tlast) := nafill(get(tlast), type = "locf"), by = id_cols]
  dt[, (var) := nafill(get(var), type = "locf"), by = id_cols]
  dt[, (var) := fifelse(
    is.na(get(tlast)) | (get(time_col) - get(tlast)) > cap_hr,
    as.numeric(NA), as.numeric(get(var))
  )]
  dt[, (tlast) := NULL]
  invisible(NULL)
}

impute_rate_start0_locf <- function(dt, var, eps = 0,
                                    id_cols = c("src", "all_id"),
                                    time_col = "time_hr") {
  if (!var %in% names(dt)) return(invisible(NULL))
  suppressWarnings(dt[, (var) := as.numeric(get(var))])
  start_tbl <- dt[!is.na(get(var)) & get(var) > eps,
                  .(start_t = suppressWarnings(min(get(time_col), na.rm = TRUE))),
                  by = id_cols]
  dt[start_tbl, on = id_cols, start_t := i.start_t]
  dt[is.na(start_t) & is.na(get(var)), (var) := 0]
  dt[!is.na(start_t) & get(time_col) < start_t & is.na(get(var)), (var) := 0]
  dt[, (var) := nafill(get(var), type = "locf"), by = id_cols]
  dt[is.na(get(var)), (var) := 0]
  dt[, start_t := NULL]
  invisible(NULL)
}

# ==============================================================================
# 10. Run imputation
# ==============================================================================

dt <- as.data.table(ts_all_ts)
dt[, src := as.character(src)]
dt[, all_id := as.character(all_id)]
dt[, time_hr := as.numeric(time_hr)]
setorderv(dt, c("src", "all_id", "time_hr"))

for (v in flag_vars) dt[, (v) := as.logical(get(v))]

for (v in static_num) {
  suppressWarnings(dt[, (v) := as.numeric(get(v))])
  dt[, (v) := nafill(get(v), type = "locf"), by = id_cols]
  dt[, (v) := nafill(get(v), type = "nocb"), by = id_cols]
}

if ("sex" %in% names(dt)) {
  sex_tbl <- dt[!is.na(sex) & sex != "", .(sex = sex[1]), by = id_cols]
  dt[sex_tbl, on = id_cols, sex := i.sex]
  dt[is.na(sex) | sex == "", sex := "Unknown"]
}

for (v in rate_vars) impute_rate_start0_locf(dt, v, eps = PARAM$eps_rate)
for (v in vitals_vars) impute_locf_cap_keepna(dt, v, cap_hr = PARAM$cap_vitals_hr)
for (v in fastlab_vars) impute_locf_cap_keepna(dt, v, cap_hr = PARAM$cap_fastlab_hr)
for (v in slowlab_vars) impute_locf_cap_keepna(dt, v, cap_hr = PARAM$cap_slowlab_hr)
for (v in sofa_vars) impute_locf_cap_keepna(dt, v, cap_hr = PARAM$cap_sofa_hr)

ts_imputed <- as_tibble(dt)

# write_rds(ts_imputed, file.path(DATA_DIR, "septic_shock_imputed_df_new"))

cat("Data extraction and imputation complete.\n")
cat("Rows:", nrow(ts_imputed), "| Patients:",
    n_distinct(paste(ts_imputed$src, ts_imputed$all_id)), "\n")
