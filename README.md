
<!-- README.md is generated from README.Rmd. Please edit that file -->

# -colorectal-surrogacy-analysis

<!-- badges: start -->
<!-- badges: end -->

This repository contains all R-code and output for the analysis of the
colorectal cancer data. The aim of this analysis is to illustrate the
evaluation of time-to-event surrogates for time-to-event true endpoints
in the causal-inference framework.

# R-code

To reproduce the results reported in Stijven et al. (2023), first run
`colorectal_model_fit.R`. The code in this file fits a set of vine
copula models, selects the best fitting model based in the Akaike
information criterion, and conducts the Monte Carlo sensitivity analysis
for the individual causal association. The results of this sensitivity
analysis are then saved to `sensitivity_analysis_results.RData`. These
results are further analyzed (plots, tables, etc.) in
`colorectal_processing.R`. Plots are saved to the `Figures` directory.

The `Code - Combining Conditions` directory contains R-code for a
supplementary analysis where a more fine-grained use of unverifiable
assumptions is illustrated. This supplemnentary analysis is described in
detail in Stijven et al. (2023, Supplementary Materials).

# Data

The `COLO_DAT.csv` file contains the colorectal cancer data. These are
pooled data from two randomized trials in colorectal cancer. In the
first trial, treatment with fluorouracil (5FU) plus interferon and
treatment with 5FU plus folinic acid were compared (Corfu-A Study Group,
1995). In the second trial, treatment with 5FU plus interferon and
treatment with 5FU alone were compared (Greco and others, 1996).

# References

Stijven, F., Alonso, A., Molenberghs, G., Van Der Elst, W., Van
Keilegom, I. (2023). An information-theoretic approach to the evaluation
of time-to-event surrogates for time-to-event true endpoints based on
causal inference.

Corfu-A Study Group. (1995). Phase iii randomized study of two
fluorouracil combinations with either interferon alfa-2a or leucovorin
in advanced colorectral cancer. J Clin Oncol 12, 921–928.

Greco, FA, Figlin, R, York, M, Einhorn, L, Schilsky, R, Marshall, EM,
Buys, SS, Froimtchuk, MJ, Schuller, J, Schuchter, L and others. (1996).
Phase iii randomized study to compare interferon alfa-2a in combination
with fluorouracil versus fluorouracil alone in patients with advanced
colorectal cancer. Journal of clinical oncology 14(10), 2674–2681.
