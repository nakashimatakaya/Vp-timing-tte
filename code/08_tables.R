################################################################################
# 08_tables.R
#
# Generate tables for the manuscript:
#   Table 1:   Baseline characteristics by observed vasopressin timing
#   S-Table 5: Clone-expanded cohort before/after steroid-adjusted IPCW weighting
#   S-Table 7: Catecholamine exposure (steroid-adjusted weighted)
#
# Requires: pt_vp6, ccw_with_outcomes, ts_out, cap99_s
# Outputs:  Word documents with formatted tables
################################################################################

library(tidyverse)
library(gtsummary)
library(flextable)
library(officer)
library(survey)

DAY28_HR <- 28 * 24

DOCX_TABLE1 <- "/Users/nakashimatakaya/Desktop/SID/Outcome/VP_steroid_table1.docx"
DOCX_STBL5  <- "/Users/nakashimatakaya/Desktop/SID/Outcome/VP_steroid_stable5.docx"
DOCX_STBL7  <- "/Users/nakashimatakaya/Desktop/SID/Outcome/VP_steroid_stable7.docx"

# ==============================================================================
# 1. Table 1: Baseline characteristics
# ==============================================================================

cat("=== Table 1 ===\n")

pt_vp6 <- pt_vp6 %>%
  mutate(
    steroid_pre_f = factor(
      as.integer(coalesce(as.logical(steroid_pre), FALSE)),
      levels = c(0L, 1L), labels = c("No", "Yes")),
    obs_group = case_when(
      tvp_rel <= 3 ~ "Ultra-early (0\u20133 h)",
      tvp_rel <= 6 ~ "Early (>3\u20136 h)",
      TRUE ~ NA_character_),
    mech_vent0_bin = factor(
      as.integer(as.character(mech_vent0) == "invasive"),
      levels = c(0L, 1L), labels = c("No", "Yes"))
  )

tbl1_data <- pt_vp6 %>%
  filter(!is.na(obs_group)) %>%
  select(obs_group, age, sex, bmi, wt0, map0, lact0, sofa0, crea0,
         norepi_equiv0, mech_vent0_bin, steroid_pre_f,
         any_of(c("death28", "aki2_28", "aki3_28", "rrt_28", "crrt_28",
                  "antiarr_28", "shockres_28", "netneg_28")))

label_list <- list(
  age = "Age, years", sex = "Sex",
  bmi = "Body mass index, kg/m\u00b2", wt0 = "Body weight, kg",
  map0 = "Mean arterial pressure, mmHg", lact0 = "Lactate, mmol/L",
  sofa0 = "SOFA score", crea0 = "Creatinine, mg/dL",
  norepi_equiv0 = "Norepinephrine-equivalent dose, \u00b5g/kg/min",
  mech_vent0_bin = "Invasive mechanical ventilation",
  steroid_pre_f = "Systemic corticosteroid use before time zero",
  death28 = "Death by day 28",
  aki2_28 = "AKI (KDIGO stage \u22652) by day 28",
  aki3_28 = "AKI (KDIGO stage \u22653) by day 28",
  rrt_28 = "RRT initiation by day 28",
  crrt_28 = "CRRT initiation by day 28",
  antiarr_28 = "Medically treated arrhythmia by day 28",
  shockres_28 = "Shock resolution by day 28",
  netneg_28 = "Net negative fluid balance by day 28"
)

tbl1 <- tbl1_data %>%
  tbl_summary(
    by = obs_group,
    label = label_list[intersect(names(label_list), names(tbl1_data))],
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
    missing = "no"
  ) %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>%
  bold_labels()

ft1 <- as_flex_table(tbl1) %>% fontsize(size = 9, part = "all") %>% autofit()

doc1 <- read_docx() %>%
  body_add_par("Table 1. Baseline characteristics and 28-day outcomes by observed vasopressin initiation timing",
               style = "heading 1") %>%
  body_add_par(paste0("Continuous variables are presented as median (interquartile range); ",
                      "categorical variables as n (%). ",
                      "28-day outcomes are crude (unweighted) event counts."),
               style = "Normal") %>%
  body_add_flextable(ft1)
print(doc1, target = DOCX_TABLE1)
cat("Table 1 saved to:", DOCX_TABLE1, "\n")

# ==============================================================================
# 2. S-Table 5: Weighted baseline characteristics
# ==============================================================================

cat("\n=== S-Table 5 ===\n")

stbl5_data <- ccw_with_outcomes %>%
  mutate(w_trunc = pmax(1e-6, pmin(weight, cap99_s)),
         arm_label = if_else(arm == "ultra", "Ultra-early (0\u20133 h)", "Early (>3\u20136 h)"))

stbl5_needed <- intersect(
  c("age", "sex", "bmi", "wt0", "map0", "lact0", "sofa0", "crea0",
    "norepi_equiv0", "mech_vent0", "steroid_pre"),
  names(pt_vp6))
stbl5_missing <- setdiff(stbl5_needed, names(stbl5_data))
if (length(stbl5_missing) > 0) {
  stbl5_data <- stbl5_data %>%
    left_join(pt_vp6 %>% select(uid, all_of(stbl5_missing)) %>% distinct(), by = "uid")
}

stbl5_data <- stbl5_data %>%
  mutate(
    steroid_pre = factor(as.integer(coalesce(as.logical(steroid_pre), FALSE)),
                         levels = c(0L, 1L), labels = c("No", "Yes")),
    mech_vent0_bin = factor(as.integer(as.character(mech_vent0) == "invasive"),
                            levels = c(0L, 1L), labels = c("No", "Yes")))

stbl5_vars <- intersect(
  c("age", "sex", "bmi", "wt0", "map0", "lact0", "sofa0", "crea0",
    "norepi_equiv0", "mech_vent0_bin", "steroid_pre"),
  names(stbl5_data))

stbl5_labels <- list(
  age = "Age, years", sex = "Sex",
  bmi = "Body mass index, kg/m\u00b2", wt0 = "Body weight, kg",
  map0 = "Mean arterial pressure, mmHg", lact0 = "Lactate, mmol/L",
  sofa0 = "SOFA score", crea0 = "Creatinine, mg/dL",
  norepi_equiv0 = "Norepinephrine-equivalent dose, \u00b5g/kg/min",
  mech_vent0_bin = "Invasive mechanical ventilation",
  steroid_pre = "Systemic corticosteroid use before time zero")
stbl5_labels <- stbl5_labels[intersect(names(stbl5_labels), stbl5_vars)]

# Unweighted
tbl_unwt <- stbl5_data %>%
  select(arm_label, all_of(stbl5_vars)) %>%
  tbl_summary(
    by = arm_label, label = stbl5_labels,
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
    missing = "no"
  ) %>%
  modify_header(label = "**Variable**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Unweighted**") %>%
  bold_labels()

# Weighted
svy_data <- stbl5_data %>%
  select(arm_label, w_trunc, all_of(stbl5_vars)) %>%
  filter(!is.na(arm_label))
svy_design <- svydesign(ids = ~1, weights = ~w_trunc, data = svy_data)

tbl_wt <- svy_design %>%
  tbl_svysummary(
    by = arm_label, include = all_of(stbl5_vars), label = stbl5_labels,
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
    missing = "no"
  ) %>%
  modify_header(label = "**Variable**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Weighted (steroid-adjusted IPCW)**") %>%
  bold_labels()

tbl_merged <- tbl_merge(
  tbls = list(tbl_unwt, tbl_wt),
  tab_spanner = c("**Unweighted**", "**Weighted (steroid-adjusted IPCW)**"))

ft_merged <- as_flex_table(tbl_merged) %>% fontsize(size = 9, part = "all") %>% autofit()

doc5 <- read_docx() %>%
  body_add_par("Supplementary Table 5. Baseline characteristics of the clone-expanded cohort before and after steroid-adjusted IPCW weighting",
               style = "heading 1") %>%
  body_add_par(paste0("Continuous variables are presented as median (interquartile range); ",
                      "categorical variables as n (%). ",
                      "Weights were estimated from an IPCW model that included systemic corticosteroid use ",
                      "before time zero as a baseline covariate."),
               style = "Normal") %>%
  body_add_flextable(ft_merged)
print(doc5, target = DOCX_STBL5)
cat("S-Table 5 saved to:", DOCX_STBL5, "\n")

# ==============================================================================
# 3. S-Table 7: Catecholamine AUC + VFD
# ==============================================================================

cat("\n=== S-Table 7 ===\n")

auc_vars <- intersect(c("ne_auc_24h", "ne_max_24h", "ne_auc_72h", "ne_max_72h"),
                      names(ts_out))

if (length(auc_vars) > 0) {
  stbl7_data <- ccw_with_outcomes %>%
    mutate(w = pmax(1e-6, pmin(weight, cap99_s)),
           arm_label = if_else(arm == "ultra", "Ultra-early (0\u20133 h)", "Early (>3\u20136 h)")) %>%
    left_join(ts_out %>% select(uid, all_of(auc_vars)) %>% distinct(), by = "uid")

  stbl7_labels <- list(
    ne_auc_24h = "NE-equivalent AUC 0\u201324 h, \u00b5g/kg",
    ne_max_24h = "NE-equivalent max 0\u201324 h, \u00b5g/kg/min",
    ne_auc_72h = "NE-equivalent AUC 0\u201372 h, \u00b5g/kg",
    ne_max_72h = "NE-equivalent max 0\u201372 h, \u00b5g/kg/min")
  stbl7_labels <- stbl7_labels[intersect(names(stbl7_labels), auc_vars)]

  svy7 <- svydesign(ids = ~1, weights = ~w,
                    data = stbl7_data %>% select(arm_label, w, all_of(auc_vars)))

  tbl7 <- svy7 %>%
    tbl_svysummary(
      by = arm_label, include = all_of(auc_vars), label = stbl7_labels,
      statistic = all_continuous() ~ "{median} ({p25}, {p75})",
      digits = all_continuous() ~ 1, missing = "no"
    ) %>%
    modify_header(label = "**Variable**") %>%
    bold_labels()

  ft7 <- as_flex_table(tbl7) %>% fontsize(size = 9, part = "all") %>% autofit()

  doc7 <- read_docx() %>%
    body_add_par("Supplementary Table 7. Catecholamine exposure (steroid-adjusted weighted)",
                 style = "heading 1") %>%
    body_add_flextable(ft7)
  print(doc7, target = DOCX_STBL7)
  cat("S-Table 7 saved to:", DOCX_STBL7, "\n")
} else {
  cat("No AUC variables found, skipping S-Table 7.\n")
}

cat("\n===== All tables complete =====\n")
