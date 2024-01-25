
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Data Application: Advanced Colorectal Cancer

<!-- badges: start -->
<!-- badges: end -->

This repository contains all R-code and output for the analysis of the
advanced colorectal cancer data. The aim of this analysis is to
illustrate the evaluation of a (composite) time-to-event surrogate for a
time-to-event true endpoint in the causal-inference framework.

The surrogacy evaluation methods are implemented in the Surrogate R-package 
on CRAN. However, the most recent version of the package can be installed
from GitHub https://github.com/florianstijven/Surrogate-development. It is
advised to install directly from GitHub for the most recent version of the 
package.

# R-code

To reproduce the results reported in Stijven et al. (2023), first run
`colorectal-model-fit-main.R`. The code in this file fits a set of vine
copula models, selects the best fitting model based on the Akaike
information criterion, and conducts the sensitivity analysis for the
individual causal association (ICA). The results of this sensitivity
analysis are then saved to `sensitivity-analysis-results-main.csv`.
These results are further analyzed (plots, tables, etc.) in
`colorectal-main-processing.R`. Plots are saved to the `Figures/`
directory.

The `VSC/R scripts/` directory contains R scripts for additional
analyses under different assumptions than the main analyses. The
corresponding results are reported in the Supplementary Materials of
Stijven et al. (2023). The results of these additional analyses are
saved in .csv-files in `VSC/results/` and further analyzed in the R
scripts in `additional-analyses-processing/`. The scripts in
`VSC/R scripts/` were run on HPC infrastructure. This HPC infrastructure
uses the slurm job scheduler. The corresponding (slurm) job files are
saved in `VSC/job-files/`. Note that these additional analyses can still
be run on a recent laptop. One sensitivity analysis takes about 90
minutes on a 10 core laptop. To run these R scripts, the `.libPaths()`
call should be commented out or changed to represent the path to the
user’s R library.

# Data

The `COLO_DAT.csv` file contains the advanced colorectal cancer data.
These are pooled data from two randomized trials in colorectal cancer.
In the first trial, treatment with fluorouracil (5FU) plus interferon
and treatment with 5FU plus folinic acid were compared (Corfu-A Study
Group, 1995). In the second trial, treatment with 5FU plus interferon
and treatment with 5FU alone were compared (Greco and others, 1996).

# References

Stijven, F., Molenberghs, G., Van Keilegom, I., Van Der Elst, W.,
Alonso, A. (2023). An information-theoretic approach to the evaluation
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
