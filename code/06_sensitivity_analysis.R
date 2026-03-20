################################################################################
# 06_sensitivity_analysis.R
#
# Sensitivity analyses:
# - 56 alternative time-zero definitions (NE threshold x sustainment rule x mode)
# - Weight truncation at 95th percentile
# - Outcome-regression-augmented model
#
# Requires: sens_results (pre-computed), ccw_tbl_steroid, cap99_s, cap95_s
# Note: sens_results was generated in a separate long-running script.
#       This file loads the cached results and prepares them for plotting.
################################################################################

library(tidyverse)
library(survival)
library(readr)

# ==============================================================================
# 1. Load pre-computed sensitivity results
# ==============================================================================

sens_results <- readr::read_rds(
 "/Users/nakashimatakaya/Desktop/SID/Data/Analysis_Data/cache_vp_secondary_20260223/sens_results_selfcontained.rds"
)

cat("Sensitivity results loaded:", nrow(sens_results), "rows\n")
cat("Columns:", paste(names(sens_results), collapse = ", "), "\n")

# ==============================================================================
# 2. Filter to cap99 and annotate
# ==============================================================================

sens_cap99 <- sens_results %>%
  filter(str_detect(model, "cap99") | cap == "cap99")

sens_cap99 <- sens_cap99 %>%
  mutate(
    facet_label = if_else(MODE == "rate", "NE infusion rate", "NE-equivalent dose"),
    y_label = sprintf("\u2265%s, K=%d, W=%d",
                      sub("^0", "", as.character(NE_THR)), SUSTAIN_K, SUSTAIN_WINDOW),
    is_main = (NE_THR == 0.25 & SUSTAIN_K == 2 & SUSTAIN_WINDOW == 2 & MODE == "equiv"),
    highlight = if_else(is_main, "Main analysis", "Alternative")
  )

cat("Main analysis HR:\n")
print(sens_cap99 %>% filter(is_main) %>% select(NE_THR, SUSTAIN_K, SUSTAIN_WINDOW, MODE, HR, lower, upper))

cat("\nRange of HRs across 56 definitions:\n")
cat(sprintf("  Min: %.3f, Max: %.3f\n", min(sens_cap99$HR), max(sens_cap99$HR)))
cat(sprintf("  All < 1.0: %s\n", all(sens_cap99$HR < 1.0)))
