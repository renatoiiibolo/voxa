# VOxA — Voxel-Aware Oxygen Model

Code for **"Voxel-aware oxygen kinetics resolves radiation-induced DNA damage retention across LET–oxygen conditions"** (Bolo & Bagunu, 2026, submitted to *Physics in Medicine and Biology*).

## Overview

VOxA predicts the oxygen enhancement ratio (OER) for DNA double-strand break (DSB) retention at particle-therapy-relevant LET and oxygen tension. It has two tiers:

- **Oxygen Model (OM):** Population-level OER via dual sigmoidal LET transitions with Z-ordering constraints and Michaelis–Menten oxygen kinetics. Calibrated on 233 OER observations from 29 published sources spanning 10 particle types.
- **Voxel-Aware (VA) extension:** Per-DSB retention probability using local energy heterogeneity from TOPAS-nBio 4.0 track-structure simulations.

## Repository structure

```
├── README.md
├── .gitignore
├── session_info.txt
│
├── step1_extract_furusawa_complete.R
├── step2_compile_literature.R
├── step3_particle_specific_dual_transitions.R
├── step4_comprehensive_validation.R
├── step5_comprehensive_validation_extended.R
├── step6_uncertainty_analysis.R
├── step7_bootstrap_ci_analysis.R
├── step8_cross_validation.R
├── step9_recalibrate_voxel_aware.R
├── step10_voxa_va_diagnostics.R
├── step11_voxa_scaling_validation.R
├── step12_dsb_retention_table.R
│
├── extract_energy_features_voxa.py
├── voxa_model.py
└── voxa_manuscript_figures.py
```

## Running the analysis

Run R scripts in order from the project root:

```bash
Rscript step1_extract_furusawa_complete.R
Rscript step2_compile_literature.R
# ... through step12
```

Before running the manuscript figure script, export model parameters from R (see the docstring at the top of `voxa_manuscript_figures.py`). Then:

```bash
python voxa_manuscript_figures.py
```

Step 7 (bootstrap) takes 20–60 minutes. All other steps are fast.

TOPAS-nBio energy features for the VA extension are pre-extracted in `voxa_features_output_calibration/` and `voxa_features_output_validation/`. To regenerate them from raw simulation output, use `extract_energy_features_voxa.py`.

## Requirements

**R** (≥ 4.3):

```r
install.packages(c("tidyverse", "MASS", "minpack.lm", "boot",
                   "jsonlite", "gridExtra", "car", "cowplot",
                   "scales", "viridis"))
```

**Python** 3.14:

```
numpy
pandas
scipy
matplotlib
```

See `session_info.txt` for exact package versions.

## Data sources

The 233-observation calibration dataset is compiled from the sources listed in Appendix C of the manuscript. Furusawa et al. (2000) data (He, C, Ne ions; V79 and HSG cells) were extracted directly from published D₁₀ tables. Tinganelli et al. (2015) data for intermediate oxygen levels are in `calibration_data/`.

## Contact

Renato III Fernan Bolo  
Department of Physical Sciences and Mathematics, University of the Philippines Manila  
rfbolo@up.edu.ph
