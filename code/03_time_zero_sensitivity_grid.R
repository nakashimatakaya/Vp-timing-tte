suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(officer)
  library(rvg)
  library(grid)
})

ROOT <- normalizePath(getwd())
SENS_RDS <- Sys.getenv(
  "VP_TIMEZERO_SENS_RDS",
  unset = file.path("path", "to", "sens_results_selfcontained.rds")
)
ASSETS <- Sys.getenv("VP_OUTPUT_DIR", unset = file.path(ROOT, "output", "time_zero_sensitivity"))
dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Times New Roman"
COL_BLACK <- "#111111"
COL_GRID <- "#E3E3E3"
COL_HIGHLIGHT <- "#FFF2A8"

fmt_hr <- function(hr, lo, hi) sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)

sens_raw <- readr::read_rds(SENS_RDS)

sens_plot <- sens_raw %>%
  filter(model == "cap99") %>%
  mutate(
    Variable = recode(
      MODE,
      "equiv" = "Norepinephrine-equivalent dose",
      "rate" = "Norepinephrine infusion rate alone"
    ),
    spec_label = sprintf("\u2265%.2f, K=%d, W=%d", NE_THR, SUSTAIN_K, SUSTAIN_WINDOW),
    is_primary = MODE == "equiv" & NE_THR == 0.25 & SUSTAIN_K == 2 & SUSTAIN_WINDOW == 2,
    marker = if_else(is_primary, "Primary specification", "Alternative specification"),
    HR_95CI = fmt_hr(HR, lower, upper)
  )

spec_levels <- sens_plot %>%
  distinct(NE_THR, SUSTAIN_K, SUSTAIN_WINDOW, spec_label) %>%
  arrange(NE_THR, SUSTAIN_K, SUSTAIN_WINDOW) %>%
  pull(spec_label)

sens_plot <- sens_plot %>%
  mutate(
    spec_label = factor(spec_label, levels = rev(spec_levels)),
    Variable = factor(
      Variable,
      levels = c("Norepinephrine-equivalent dose", "Norepinephrine infusion rate alone")
    )
  )

or11_timezero_plot <- ggplot(sens_plot, aes(x = HR, y = spec_label)) +
  geom_rect(
    data = sens_plot %>% filter(is_primary),
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = as.numeric(spec_label) - 0.46,
      ymax = as.numeric(spec_label) + 0.46
    ),
    inherit.aes = FALSE,
    fill = COL_HIGHLIGHT,
    color = NA
  ) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.55, color = "#555555") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18, linewidth = 0.38, color = COL_BLACK) +
  geom_point(aes(shape = marker), size = 2.2, color = COL_BLACK, fill = COL_BLACK) +
  scale_shape_manual(values = c("Alternative specification" = 16, "Primary specification" = 18)) +
  facet_wrap(~Variable, nrow = 1) +
  scale_x_continuous(breaks = seq(0.3, 1.1, 0.1), limits = c(0.30, 1.10)) +
  labs(
    x = "Hazard ratio (ultra-early vs early)",
    y = "Threshold, K, W",
    shape = NULL
  ) +
  theme_bw(base_family = FONT_FAMILY, base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.35),
    panel.grid.major.x = element_line(color = COL_GRID, linewidth = 0.35),
    panel.border = element_rect(color = COL_BLACK, linewidth = 0.45),
    strip.background = element_rect(fill = "white", color = COL_BLACK, linewidth = 0.45),
    strip.text = element_text(face = "bold", size = 12.5),
    axis.text.y = element_text(size = 7.7, color = COL_BLACK),
    axis.text.x = element_text(size = 8.8, color = COL_BLACK),
    axis.title.x = element_text(size = 12.5, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 11.5, face = "bold", margin = margin(r = 7)),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 10.5),
    legend.margin = margin(b = 1),
    plot.margin = margin(6, 8, 6, 6)
  )

readr::write_csv(
  sens_plot %>%
    arrange(Variable, NE_THR, SUSTAIN_K, SUSTAIN_WINDOW) %>%
    transmute(
      variable = as.character(Variable),
      threshold = NE_THR,
      K = SUSTAIN_K,
      W = SUSTAIN_WINDOW,
      HR,
      lower,
      upper,
      n_uid,
      events,
      HR_95CI,
      primary_specification = is_primary
    ),
  file.path(ASSETS, "OR11_sensitivity_estimates.csv")
)
readr::write_csv(
  sens_plot %>%
    arrange(Variable, NE_THR, SUSTAIN_K, SUSTAIN_WINDOW),
  file.path(ASSETS, "OR11_timezero_sensitivity_56_estimates.csv")
)

ggsave(
  file.path(ASSETS, "OR11_sensitivity_forest.pdf"),
  or11_timezero_plot,
  width = 11.2,
  height = 6.2,
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(ASSETS, "OR11_sensitivity_forest.png"),
  or11_timezero_plot,
  width = 11.2,
  height = 6.2,
  dpi = 320,
  bg = "white"
)

pptx_path <- file.path(ASSETS, "OR11_timezero_sensitivity_editable.pptx")
ppt <- read_pptx()
ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
ppt <- ph_with(
  ppt,
  "Online Resource 11. Sensitivity of the primary outcome to alternative time-zero definitions",
  location = ph_location(left = 0.25, top = 0.10, width = 9.5, height = 0.35),
  gp = fp_text(font.family = FONT_FAMILY, font.size = 12, bold = TRUE)
)
ppt <- ph_with(
  ppt,
  dml(ggobj = or11_timezero_plot),
  location = ph_location(left = 0.20, top = 0.52, width = 9.60, height = 5.35)
)
print(ppt, target = pptx_path)

deck_path <- file.path(ASSETS, "VP_online_resources_7_16_figures_editable.pptx")
if (file.exists(deck_path)) {
  deck <- read_pptx(deck_path)
  deck_size <- slide_size(deck)
  deck <- on_slide(deck, index = 4)
  deck <- ph_with(
    deck,
    dml(code = grid.rect(gp = gpar(fill = "white", col = "white"))),
    location = ph_location(left = 0, top = 0, width = deck_size$width, height = deck_size$height)
  )
  deck <- ph_with(
    deck,
    "Online Resource 11. Sensitivity of the primary outcome to alternative time-zero definitions",
    location = ph_location(left = 0.25, top = 0.10, width = 9.5, height = 0.35),
    gp = fp_text(font.family = FONT_FAMILY, font.size = 12, bold = TRUE)
  )
  deck <- ph_with(
    deck,
    dml(ggobj = or11_timezero_plot),
    location = ph_location(left = 0.20, top = 0.52, width = 9.60, height = 5.35)
  )
  print(deck, target = deck_path)
}

message("Wrote corrected OR11 assets to: ", ASSETS)
