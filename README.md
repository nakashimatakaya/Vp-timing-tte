# Ultra-early versus early adjunctive vasopressin initiation after norepinephrine escalation in septic shock

This repository contains the analysis code used for the revised target trial emulation comparing ultra-early (0-3 hours) versus early (>3-6 hours) adjunctive vasopressin initiation after norepinephrine escalation in patients with septic shock.

Only analysis-stage code is provided here. Database extraction, harmonization, and preprocessing scripts are intentionally not included, because the source datasets require credentialed access through PhysioNet and local preprocessing produces derived RDS files that cannot be redistributed.

## Data sources

The study used two publicly available critical care databases:

- Medical Information Mart for Intensive Care IV (MIMIC-IV), version 3.1
- eICU Collaborative Research Database (eICU-CRD), version 2.0

No patient-level data or derived analysis datasets are included in this repository.

## Repository structure

```text
vp-timing-tte/
├── README.md
├── LICENSE
├── .gitignore
└── code/
    ├── 01_primary_revision_analysis.Rmd
    ├── 02_final_tables_figures.R
    └── 03_time_zero_sensitivity_grid.R
```

## Code contents

- `01_primary_revision_analysis.Rmd`  
  Main revision analysis, including clone-censor-weight construction, inverse-probability-of-censoring-weighted Cox models, baseline-only IPCW sensitivity analysis, norepinephrine-restricted analyses, inotrope-related sensitivity analyses, secondary outcomes, subgroup analyses, descriptive tables, risk difference, restricted mean survival time, and response-letter summaries.

- `02_final_tables_figures.R`  
  Final manuscript and supplementary table/figure generation based on the revised norepinephrine time-zero analysis, including the main figures, primary and key sensitivity analyses, subgroup forest plot, Table 1, and supplementary descriptive tables.

- `03_time_zero_sensitivity_grid.R`  
  Forest plot and table generation for the grid of alternative time-zero definitions using 56 combinations of dose threshold, sustainment duration, and verification window.

## Required local inputs

The scripts expect preprocessed analysis RDS files produced from MIMIC-IV and eICU-CRD. These files are not included. Set the following environment variables before running the scripts:

```r
Sys.setenv(
  VP_ANALYSIS_DATA_DIR = "/path/to/analysis_data",
  VP_ELIGIBLE_PATIENT_RDS = "/path/to/pt_all_eligible.rds",
  VP_TS_IMPUTED_RDS = "/path/to/septic_shock_imputed_df_new",
  VP_FLUID_HOURLY_RDS = "/path/to/fluid_hourly.rds",
  VP_TIMEZERO_SENS_RDS = "/path/to/sens_results_selfcontained.rds",
  VP_OUTPUT_DIR = "/path/to/output"
)
```

The expected derived files include, among others, the `revICM` analysis objects used by the revision analyses:

- `vp_ts_out_revICM_20260517.rds`
- `comorbidity_tbl.rds`
- `safety_tbl.rds`
- `ne_flags_tbl.rds`
- `t0_neonly_tbl.rds`
- `inotrope_flags_tbl.rds`
- `vp_dose_tbl.rds`
- `vp_interrupt_tbl.rds`
- `fluid_summary_tbl.rds`
- `steroid_timing_tbl.rds`

## Requirements

R version 4.3.2 or later is recommended.

Core R packages:

```r
install.packages(c(
  "tidyverse", "survival", "scales", "patchwork", "splines",
  "gtsummary", "flextable", "officer", "rvg", "grid",
  "sandwich", "fs"
))
```

## Reproduction notes

Run the scripts after the local preprocessed RDS files have been created and the environment variables above have been set. The scripts are intended to document and reproduce the analysis stage, not to recreate the credentialed database extraction process.

## Citation

Nakashima T, Ichinomiya T, Nakajima M, Shinozaki T, Shibata J, Goto T, Sato S, Hara T. Ultra-early versus early adjunctive vasopressin initiation after norepinephrine escalation in septic shock: a target trial emulation. 2026.

## License

This project is licensed under the MIT License. See `LICENSE` for details.
