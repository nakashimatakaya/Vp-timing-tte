################################################################################
# 07_figures.R
#
# Generate all figures for the manuscript:
#   Figure 1: Study flow diagram (Gmisc/grid)
#   Figure 2: Weighted survival curves + NAR + 95% CI bands
#   Figure 3: Panel A (primary robustness) + Panel B (secondary outcomes)
#   S-Figure 3: Weight distribution
#   S-Figure 4: Love plot (SMD before/after weighting)
#   S-Figure 5: Subgroup forest plot
#   S-Figure 6: Database-specific estimates
#   S-Figure 7: 56 time-zero sensitivity (faceted, highlighted)
#   S-Figure 8: Weighted NE trajectory 0-72h
#
# Requires: Objects from scripts 03-06
# Outputs:  PowerPoint file with all figures
################################################################################

library(tidyverse)
library(survival)
library(scales)
library(patchwork)
library(officer)
library(rvg)
library(Gmisc)
library(grid)

DAY28_HR  <- 28 * 24
PPTX_PATH <- "/Users/nakashimatakaya/Desktop/SID/Outcome/VP_steroid_adjusted_figures.pptx"

fmt_n <- function(x) format(as.integer(x), big.mark = ",")

# ==============================================================================
# Figure 1: Study flow diagram
# ==============================================================================

plot_fig1 <- function() {
  grid.newpage()
  txt_gp  <- gpar(fontsize = 14)
  box_gp  <- gpar(fill = "white")

  box1 <- boxGrob(
    paste0("Adults with septic shock in MIMIC-IV and eICU-CRD\n",
           "(N = 21,587; MIMIC-IV 14,531, eICU-CRD 7,056)"),
    x = 0.42, y = 0.85, txt_gp = txt_gp, box_gp = box_gp)
  box_excl <- boxGrob(
    paste0("Excluded (n = 16,430)\n",
           "NE-equivalent dose <0.25 \u00b5g/kg/min\n",
           "Vasopressin before time zero\n",
           "No post-escalation follow-up"),
    x = 0.82, y = 0.735, txt_gp = txt_gp, box_gp = box_gp)
  box2 <- boxGrob(
    paste0("Eligible at time zero\n",
           "NE-equivalent dose \u22650.25 \u00b5g/kg/min, no prior vasopressin\n",
           "(n = 5,157; MIMIC-IV 2,721, eICU-CRD 2,436)"),
    x = 0.42, y = 0.61, txt_gp = txt_gp, box_gp = box_gp)
  box3 <- boxGrob(
    paste0("Initiated vasopressin within 6 h after time zero\n",
           "(n = 1,266; MIMIC-IV 624, eICU-CRD 642)"),
    x = 0.42, y = 0.45, txt_gp = txt_gp, box_gp = box_gp)
  box_u <- boxGrob(
    paste0("Ultra-early (0\u20133 h)\n",
           "(n = 943; MIMIC-IV 422, eICU-CRD 521)"),
    x = 0.22, y = 0.25, txt_gp = txt_gp, box_gp = box_gp)
  box_e <- boxGrob(
    paste0("Early (>3\u20136 h)\n",
           "(n = 323; MIMIC-IV 202, eICU-CRD 121)"),
    x = 0.62, y = 0.25, txt_gp = txt_gp, box_gp = box_gp)

  connectGrob(box1, box2, type = "vertical")
  connectGrob(box1, box_excl, type = "L")
  connectGrob(box2, box3, type = "vertical")
  connectGrob(box3, box_u, type = "N")
  connectGrob(box3, box_e, type = "N")

  box1; box_excl; box2; box3; box_u; box_e
}

# ==============================================================================
# Figure 2: Weighted survival curves with CI bands and NAR
# ==============================================================================

build_fig2 <- function() {
  sf_df <- tibble(
    time = sf$time, surv = sf$surv, upper = sf$upper, lower = sf$lower,
    arm = rep(names(sf$strata), sf$strata)
  ) %>%
    mutate(
      arm = str_replace(arm, "arm=", ""),
      arm = case_when(arm == "ultra" ~ "0\u20133 h (ultra-early)",
                      arm == "early" ~ ">3\u20136 h (early)", TRUE ~ arm),
      arm = factor(arm, levels = c("0\u20133 h (ultra-early)", ">3\u20136 h (early)"))
    )

  nar_times <- seq(0, DAY28_HR, by = 168)
  nar_df <- map_dfr(c("ultra", "early"), function(a) {
    d <- tmp_surv %>% filter(arm == a)
    label <- if_else(a == "ultra", "0\u20133 h", ">3\u20136 h")
    map_dfr(nar_times, function(t_hr) {
      tibble(arm_label = label, time = t_hr,
             n_risk = sum(d$tstart <= t_hr & d$tstop > t_hr, na.rm = TRUE))
    })
  })

  annot_text <- paste0(
    "Day-28 mortality risk\n",
    sprintf("  0\u20133 h: %.1f%%\n", risk_ultra * 100),
    sprintf("  >3\u20136 h: %.1f%%\n", risk_early * 100),
    sprintf("Risk difference: %.1f pp (%.1f to %.1f)\n", rd_point, rd_ci[1], rd_ci[2]),
    sprintf("RMST difference: +%.2f days (%.2f to %.2f)\n", rmst_diff_pt, rmst_ci[1], rmst_ci[2]),
    sprintf("HR: %.2f (%.2f\u2013%.2f)", res_99$HR, res_99$lower, res_99$upper))

  p_main <- ggplot(sf_df, aes(x = time, y = surv, linetype = arm)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = arm), alpha = 0.15, colour = NA) +
    geom_step(linewidth = 0.7, colour = "black") +
    scale_fill_manual(values = c("0\u20133 h (ultra-early)" = "grey30",
                                 ">3\u20136 h (early)" = "grey60"), guide = "none") +
    scale_linetype_manual(values = c("0\u20133 h (ultra-early)" = "solid",
                                     ">3\u20136 h (early)" = "dashed")) +
    scale_x_continuous(breaks = nar_times, labels = nar_times / 24,
                       limits = c(0, DAY28_HR), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       labels = percent_format(accuracy = 1), expand = c(0.01, 0)) +
    annotate("text", x = DAY28_HR * 0.42, y = 0.18, label = annot_text,
             hjust = 0, vjust = 0, size = 3.2, lineheight = 1.2) +
    labs(x = "Days after time zero", y = "Survival probability", linetype = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = c(0.78, 0.92), legend.background = element_blank(),
          legend.key.width = unit(1.5, "cm"), panel.grid.minor = element_blank(),
          plot.margin = margin(5.5, 10, 40, 10))

  p_nar <- ggplot(nar_df, aes(x = time, y = arm_label, label = n_risk)) +
    geom_text(size = 3.2) +
    scale_x_continuous(breaks = nar_times, limits = c(0, DAY28_HR), expand = c(0.01, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(axis.text.y = element_text(size = 9, hjust = 1),
          plot.margin = margin(0, 10, 5.5, 10))

  p_main / p_nar + plot_layout(heights = c(5, 1))
}

# ==============================================================================
# Figure 3A: Primary robustness forest
# ==============================================================================

build_fig3A <- function() {
  primary_order <- c("IPCW Cox (99th percentile truncation)",
                     "IPCW Cox (95th percentile truncation)",
                     "IPCW Cox + outcome regression")
  df <- res_primary %>%
    mutate(hr_label = sprintf("%.2f (%.2f\u2013%.2f)", HR, lower, upper),
           label = factor(model, levels = rev(primary_order)))
  ggplot(df, aes(x = HR, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15, linewidth = 0.6) +
    geom_point(size = 3.5) +
    geom_text(aes(label = hr_label), vjust = -1.1, size = 3.6, colour = "grey20") +
    scale_x_log10(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
    coord_cartesian(xlim = c(0.48, 1.08), clip = "off") +
    labs(title = "A  Primary outcome", x = "Hazard ratio (ultra-early vs early)", y = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold"), plot.margin = margin(24, 10, 5.5, 5.5))
}

# ==============================================================================
# Figure 3B: Secondary outcomes forest
# ==============================================================================

build_fig3B <- function() {
  df <- sec_results %>%
    filter(!is.na(HR)) %>%
    mutate(hr_label = sprintf("%.2f (%.2f\u2013%.2f)", HR, LCL, UCL),
           label = factor(outcome, levels = rev(sec_labels)))
  x_lo <- min(0.85, min(df$LCL, na.rm = TRUE) * 0.95)
  x_hi <- max(1.10, max(df$UCL, na.rm = TRUE) * 1.05)
  ggplot(df, aes(x = HR, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.20, linewidth = 0.50) +
    geom_point(size = 3) +
    geom_text(aes(label = hr_label), vjust = -1.1, size = 3.0, colour = "grey20") +
    scale_x_log10() +
    coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
    labs(title = "B  Secondary outcomes",
         x = "Subdistribution hazard ratio (ultra-early vs early)", y = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold"), plot.margin = margin(24, 12, 5.5, 5.5))
}

# ==============================================================================
# S-Figure 3: Weight distribution
# ==============================================================================

build_sfig3 <- function() {
  df <- ccw_with_outcomes %>%
    mutate(w_trunc = pmin(weight, cap99_s),
           arm_label = if_else(arm == "ultra", "0\u20133 h (ultra-early)", ">3\u20136 h (early)"))
  ggplot(df, aes(x = w_trunc, fill = arm_label)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("0\u20133 h (ultra-early)" = "grey30", ">3\u20136 h (early)" = "grey70")) +
    labs(x = "Stabilized IPCW (truncated at 99th percentile)", y = "Count", fill = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "top")
}

# ==============================================================================
# S-Figure 4: Love plot (SMD)
# ==============================================================================

build_sfig4 <- function() {
  var_name_map <- c(
    age = "Age, years", sex = "Male sex",
    bmi = "Body mass index, kg/m\u00b2", wt0 = "Body weight, kg",
    map0 = "Mean arterial pressure, mmHg", lact0 = "Lactate, mmol/L",
    sofa0 = "SOFA score", crea0 = "Creatinine, mg/dL",
    norepi_equiv0 = "Norepinephrine-equivalent dose, \u00b5g/kg/min",
    mech_vent0 = "Invasive mechanical ventilation",
    steroid_pre = "Systemic corticosteroid use")
  cov_vars <- intersect(names(var_name_map), names(ccw_with_outcomes))
  ccw_num <- ccw_with_outcomes %>%
    mutate(across(any_of(c("sex", "mech_vent0")), ~ as.numeric(as.factor(.x))),
           w_trunc = pmax(1e-6, pmin(weight, cap99_s)))
  calc_smd <- function(dat, var, wt = NULL) {
    d <- dat %>% filter(!is.na(.data[[var]]))
    if (nrow(d) == 0) return(NA_real_)
    x <- as.numeric(d[[var]]); a <- d$arm == "ultra"
    w <- if (!is.null(wt) && wt %in% names(d)) d[[wt]] else rep(1, nrow(d))
    m1 <- weighted.mean(x[a], w[a], na.rm = TRUE)
    m0 <- weighted.mean(x[!a], w[!a], na.rm = TRUE)
    v1 <- sum(w[a] * (x[a] - m1)^2, na.rm = TRUE) / sum(w[a], na.rm = TRUE)
    v0 <- sum(w[!a] * (x[!a] - m0)^2, na.rm = TRUE) / sum(w[!a], na.rm = TRUE)
    denom <- sqrt((v1 + v0) / 2)
    if (denom < 1e-10) return(0)
    (m1 - m0) / denom
  }
  smd_df <- tibble(variable = cov_vars) %>%
    rowwise() %>%
    mutate(smd_raw = calc_smd(ccw_num, variable),
           smd_wt = calc_smd(ccw_num, variable, "w_trunc"),
           var_label = unname(var_name_map[variable])) %>%
    ungroup() %>%
    pivot_longer(c(smd_raw, smd_wt), names_to = "type", values_to = "smd") %>%
    mutate(type = if_else(type == "smd_raw", "Before weighting", "After weighting"),
           var_label = factor(var_label, levels = rev(unname(var_name_map[cov_vars]))))
  ggplot(smd_df, aes(x = abs(smd), y = var_label, colour = type, shape = type)) +
    geom_vline(xintercept = 0.1, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3) +
    scale_colour_manual(values = c("Before weighting" = "grey50", "After weighting" = "black")) +
    labs(x = "Absolute standardized mean difference", y = NULL, colour = NULL, shape = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "top")
}

# ==============================================================================
# S-Figure 5: Subgroup forest
# ==============================================================================

build_sfig5 <- function() {
  df <- sub_results %>%
    filter(!is.na(HR)) %>%
    mutate(display = paste0(subgroup, ": ", level),
           hr_label = sprintf("%.2f (%.2f\u2013%.2f)", HR, lower, upper),
           display = fct_rev(fct_inorder(display)))
  ggplot(df, aes(x = HR, y = display)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15, linewidth = 0.5) +
    geom_text(aes(label = hr_label), vjust = -1.2, size = 3.2, colour = "grey20") +
    scale_x_log10(breaks = c(0.3, 0.5, 1.0, 2.0, 3.0, 5.0)) +
    coord_cartesian(xlim = c(0.3, 5.5)) +
    labs(x = "Hazard ratio (ultra-early vs early)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
}

# ==============================================================================
# S-Figure 6: Database-specific forest
# ==============================================================================

build_sfig6 <- function() {
  df <- db_results %>%
    filter(!is.na(HR)) %>%
    mutate(hr_label = sprintf("%.2f (%.2f\u2013%.2f)", HR, lower, upper),
           database = fct_rev(fct_inorder(database)))
  ggplot(df, aes(x = HR, y = database)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_point(size = 3.5) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15, linewidth = 0.6) +
    geom_text(aes(label = hr_label), vjust = -1.2, size = 3.4, colour = "grey20") +
    scale_x_log10(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1)) +
    coord_cartesian(xlim = c(0.4, 1.15)) +
    labs(x = "Hazard ratio (ultra-early vs early)", y = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
}

# ==============================================================================
# S-Figure 7: 56 time-zero definitions (faceted)
# ==============================================================================

build_sfig7 <- function() {
  df <- sens_cap99 %>%
    arrange(facet_label, HR) %>%
    mutate(y_label = fct_inorder(y_label))
  ggplot(df, aes(x = HR, y = y_label)) +
    geom_rect(data = df %>% filter(is_main),
              aes(xmin = -Inf, xmax = Inf,
                  ymin = as.numeric(y_label) - 0.4, ymax = as.numeric(y_label) + 0.4),
              fill = "#FFFF99", alpha = 0.5, inherit.aes = FALSE) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, linewidth = 0.3) +
    geom_point(aes(shape = highlight), size = 1.8) +
    scale_shape_manual(values = c("Main analysis" = 18, "Alternative" = 16)) +
    facet_wrap(~ facet_label, scales = "free_y") +
    scale_x_continuous(breaks = seq(0.2, 1.2, 0.2)) +
    coord_cartesian(xlim = c(0.3, 1.15)) +
    labs(x = "Hazard ratio (ultra-early vs early)", y = NULL, shape = NULL) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(face = "bold", size = 11),
          axis.text.y = element_text(size = 7),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          legend.position = "top")
}

# ==============================================================================
# S-Figure 8: Weighted NE trajectory
# ==============================================================================

build_sfig8 <- function() {
  uid_weights <- ccw_with_outcomes %>%
    group_by(uid, arm) %>%
    summarise(w = first(pmax(1e-6, pmin(weight, cap99_s))), .groups = "drop")
  ne_traj <- ts_tbl_data %>%
    mutate(uid = str_c(src, all_id, sep = "::")) %>%
    inner_join(pt_vp6 %>% select(uid, t0_hr) %>% distinct(), by = "uid") %>%
    mutate(t_rel_hr = time_hr - t0_hr) %>%
    filter(t_rel_hr >= 0, t_rel_hr <= 72, !is.na(norepi_equiv)) %>%
    inner_join(uid_weights, by = "uid")
  ne_summary <- ne_traj %>%
    mutate(t_bin = floor(t_rel_hr)) %>%
    group_by(arm, t_bin) %>%
    summarise(
      ne_mean = weighted.mean(norepi_equiv, w, na.rm = TRUE),
      ne_se = sqrt(sum(w * (norepi_equiv - weighted.mean(norepi_equiv, w, na.rm = TRUE))^2) /
                     sum(w, na.rm = TRUE)) / sqrt(n()),
      .groups = "drop") %>%
    mutate(arm_label = if_else(arm == "ultra", "0\u20133 h (ultra-early)", ">3\u20136 h (early)"))
  ggplot(ne_summary, aes(x = t_bin, y = ne_mean, linetype = arm_label)) +
    geom_ribbon(aes(ymin = ne_mean - ne_se, ymax = ne_mean + ne_se, fill = arm_label),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8, colour = "black") +
    scale_linetype_manual(values = c("0\u20133 h (ultra-early)" = "solid", ">3\u20136 h (early)" = "dashed")) +
    scale_fill_manual(values = c("0\u20133 h (ultra-early)" = "grey30", ">3\u20136 h (early)" = "grey70")) +
    labs(x = "Hours after time zero",
         y = "Norepinephrine-equivalent dose (\u00b5g/kg/min)",
         linetype = NULL, fill = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "top")
}

# ==============================================================================
# Export all figures to PowerPoint
# ==============================================================================

fig2_combined <- build_fig2()
fig3A  <- build_fig3A()
fig3B  <- build_fig3B()
fig3   <- fig3A / fig3B + plot_layout(heights = c(1, 1.15))
sfig3  <- build_sfig3()
sfig4  <- build_sfig4()
sfig5  <- build_sfig5()
sfig6  <- build_sfig6()
sfig7  <- build_sfig7()
sfig8  <- build_sfig8()

add_fig <- function(ppt, fig, title) {
  ppt %>%
    add_slide(layout = "Blank") %>%
    ph_with(dml(ggobj = fig),
            location = ph_location(left = 0.25, top = 0.25, width = 9.5, height = 7.0))
}

ppt <- read_pptx()
ppt <- add_fig(ppt, fig2_combined, "Figure 2")
ppt <- add_fig(ppt, fig3, "Figure 3")
ppt <- add_fig(ppt, sfig3, "S-Figure 3")
ppt <- add_fig(ppt, sfig4, "S-Figure 4")
ppt <- add_fig(ppt, sfig5, "S-Figure 5")
ppt <- add_fig(ppt, sfig6, "S-Figure 6")
ppt <- add_fig(ppt, sfig7, "S-Figure 7")
ppt <- add_fig(ppt, sfig8, "S-Figure 8")

print(ppt, target = PPTX_PATH)
cat("All figures saved to:", PPTX_PATH, "\n")
