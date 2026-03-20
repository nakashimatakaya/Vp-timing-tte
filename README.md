# Ultra-early versus early adjunctive vasopressin initiation after norepinephrine escalation in septic shock: a target trial emulation

## Overview

This repository contains the analysis code for a target trial emulation comparing ultra-early (0--3 hours) versus early (>3--6 hours) adjunctive vasopressin initiation after norepinephrine escalation to 0.25 µg/kg/min or higher in patients with septic shock.

The study used clone-censor-weight methods with inverse-probability-of-censoring-weighted Cox models to estimate the per-protocol effect on 28-day mortality.

## Data sources

This study used two publicly available critical care databases:

- **MIMIC-IV** (version 3.1): [https://physionet.org/content/mimiciv/](https://physionet.org/content/mimiciv/)
- **eICU-CRD** (version 2.0): [https://physionet.org/content/eicu-crd/](https://physionet.org/content/eicu-crd/)

Access requires completion of the CITI training program and a signed data use agreement via PhysioNet. **No patient-level data are included in this repository.**

## Repository structure

```
vp-timing-tte/
├── README.md
├── LICENSE
├── .gitignore
├── code/
│   ├── 01_data_extraction.R        # Cohort extraction and hourly harmonization (ricu)
│   ├── 02_steroid_extraction.R     # Systemic corticosteroid exposure extraction
│   ├── 03_clone_censor_weight.R    # Clone-censor-weight framework with steroid-adjusted IPCW
│   ├── 04_primary_analysis.R       # Primary outcome: 28-day mortality (Cox, RMST, RD)
│   ├── 05_secondary_outcomes.R     # Secondary outcomes (AKI, RRT, arrhythmia, etc.)
│   ├── 06_sensitivity_analysis.R   # 56 time-zero definitions, weight truncation, outcome regression
│   ├── 07_figures.R                # Figures 1--3 and Supplementary Figures 3--8
│   └── 08_tables.R                 # Table 1 and Supplementary Tables 5, 7
└── output/
    └── .gitkeep
```

## Requirements

### R version

R >= 4.3.2

### R packages

```r
install.packages(c(
  "tidyverse", "survival", "survey", "splines",
  "gtsummary", "flextable", "scales", "patchwork",
  "officer", "rvg", "Gmisc"
))
```

Data extraction additionally requires:

```r
install.packages("ricu")
```

## Reproduction

1. Obtain access to MIMIC-IV and eICU-CRD via [PhysioNet](https://physionet.org/).
2. Run scripts in numerical order (`01_` through `08_`).
3. Paths to local data files are defined at the top of each script and must be updated to match your environment.

## Citation

> Nakashima T, Ichinomiya T, Nakajima M, Shinozaki T, Shibata J, Goto T, Sato S, Hara T. Ultra-early versus early adjunctive vasopressin initiation after norepinephrine escalation in septic shock: a target trial emulation. *[Journal]* (2026).

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Contact

Takaya Nakashima  
Department of Anesthesiology and Intensive Care Medicine, Nagasaki University  
Email: bb55322029@ms.nagasaki-u.ac.jp
# -vp-timing-tte
