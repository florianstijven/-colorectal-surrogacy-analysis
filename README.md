
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Data Application: Advanced Colorectal Cancer

<!-- badges: start -->
<!-- badges: end -->

This repository contains all R-code and output for the analysis of the
advanced colorectal cancer data. The aim of this analysis is to
illustrate the evaluation of a (composite) time-to-event surrogate for a
time-to-event true endpoint in the information-theoretic
causal-inference framework.

This project relies on the `renv` package for managing R packages and
ensuring a reproducible environment. When opening this project in
RStudio, `renv` is automatically initialized and will suggest to
install/update packages. To ensure reproducibility, it is advised to
follow this suggestion.

# R-code

## Main Analysis

To reproduce the results reported in Stijven et al. (2024), first run
`colorectal-model-fit-main.R`. The code in this file fits a set of vine
copula models and selects the best fitting model based on the Akaike
information criterion. The fitted models are saved into
`results\fitted-models.rds`. The best fitted model is saved into
`results\best-fitted-model.rds` and is summarized in
`results\best-fitted-model-summary.txt`.

Next, run `colorectal-sensitivity-analysis.R` to conduct the sensitivity
analysis for the individual causal association (ICA), based on the best
fitted model. The results of this sensitivity analysis are then saved to

- `results\sensitivity-analysis-results-main-sicc.rds` and
  `results\sensitivity-intervals-Rh-subset.rds`: results of the
  sensitivity analysis where the ICA is quantified in terms of the
  squared informational coefficient of correlation (SICC) in the subset
  of patients that progress before dying in both treatment arms.
- `results\sensitivity-analysis-results-main-sprho.rds`,
  `results\sensitivity-intervals-sprho-full.rds`, and
  `results\sensitivity-intervals-sprho-subset.rds`: results of the
  sensitivity analysis where the ICA is quantified in terms of
  Spearman’s rho in the full population and in the subset of patients
  that progress before dying in both treatment arms.

## Analysis Under Relaxed Assumptions

The sensitivity analysis under relaxed assumptions is implemented in
`colorectal-sensitivity-analysis-relaxed.R` and the results are saved to
`results\sensitivity-analysis-results-relaxed.rds`.

## Results

All tables and figures presented in Stijven et al. (2024) are saved into
the `paper-figures-tables\` directory when running (i) the scripts
mentioned above and (ii) running `colorectal-main-processing.R`.
Additional figures and tables that are not presented in Stijven et
al. (2024) are saved into the `additional-figures-tables\` directory.

# Data

The `COLO_DAT.csv` file contains the advanced colorectal cancer data.
These are pooled data from two randomized trials in colorectal cancer.
In the first trial, treatment with fluorouracil (5FU) plus interferon
and treatment with 5FU plus folinic acid were compared (Corfu-A Study
Group, 1995). In the second trial, treatment with 5FU plus interferon
and treatment with 5FU alone were compared (Greco and others, 1996).

# Additional Files

Files that were not mentioned above are related to the `renv` package,
licensing, RStudio project files, or files related to HPC
infrastructure. Specifically, `.slurm` files contain scripts for running
certain R scripts on HPC infrastructure.

# References

Stijven, F., Molenberghs, G., Van Keilegom, I., Van Der Elst, W.,
Alonso, A. (2024). Time-to-Event Surrogates for Time-to-Event True
Endpoints: An Information-Theoretic Approach Based on Causal Inference.

Corfu-A Study Group. (1995). Phase iii randomized study of two
fluorouracil combinations with either interferon alfa-2a or leucovorin
in advanced colorectral cancer. J Clin Oncol 12, 921–928.

Greco, FA, Figlin, R, York, M, Einhorn, L, Schilsky, R, Marshall, EM,
Buys, SS, Froimtchuk, MJ, Schuller, J, Schuchter, L and others. (1996).
Phase iii randomized study to compare interferon alfa-2a in combination
with fluorouracil versus fluorouracil alone in patients with advanced
colorectal cancer. Journal of clinical oncology 14(10), 2674–2681.
