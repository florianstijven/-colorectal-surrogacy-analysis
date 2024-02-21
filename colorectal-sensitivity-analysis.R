# Setup -------------------------------------------------------------------

ncores = 10
# Number of MC replications in the sensitivity analysis
n_sim = 5e3
# Number of MC samples for computing the ICA and related measures.
n_prec = 1e4
# Number of bootstrap replications for computing uncertainty intervals.
B = 5e2

library(Surrogate)
library(tidyverse)
# Load the best fitted model.
best_fitted_model = readRDS("results/best-fitted-model.rds")

# Sensitivity Analysis ----------------------------------------------------

# The sensitivity analysis is implemented in this file, but the results of the
# sensitivity analysis are saved and processed elsewhere. This is done because
# the sensitivity analysis is computer intensive and not interactive; whereas
# processing the results is not computer intensive, but interactive.

# Run the sensitivity analysis with the main assumptions once for ICA = SICC and
# once for ICA = Spearman's rho.
set.seed(1)
a = Sys.time()
sens_results_sicc = sensitivity_analysis_SurvSurv_copula(
  fitted_model = best_fitted_model,
  n_sim = n_sim,
  n_prec = n_prec,
  ncores = ncores,
  marg_association = TRUE,
  cond_ind = TRUE,
  composite = TRUE,
  degrees = 0,
  copula_family2 = "gaussian",
  lower = c(0.5, 0, 0, 0.15),
  upper = c(0.95, 0, 0, 0.8)
)
set.seed(1)
sensitivity_intervals_Rh_subset = sensitivity_intervals_Dvine(
  fitted_model = best_fitted_model,
  sens_results = sens_results_sicc,
  B = B,
  ncores = ncores
)
set.seed(1)
sensitivity_intervals_sprho_full = sensitivity_intervals_Dvine(
  fitted_model = best_fitted_model,
  sens_results = sens_results_sicc,
  measure = "sp_rho",
  B = B,
  ncores = ncores
)
print(Sys.time() - a)
set.seed(1)

# The SICC can be replaced with any measure by replacing the mutual information
# estimator with an estimator of -0.5 * log(1 - measure); next, we use measure =
# stats::cor(x, y, method = "spearman") to compute Spearman's rho.
a = Sys.time()
sens_results_sprho = sensitivity_analysis_SurvSurv_copula(
  fitted_model = best_fitted_model,
  n_sim = n_sim,
  n_prec = n_prec,
  ncores = ncores,
  marg_association = TRUE,
  cond_ind = TRUE,
  composite = TRUE,
  degrees = 0,
  copula_family2 = "gaussian",
  lower = c(0.5, 0, 0, 0.15),
  upper = c(0.95, 0, 0, 0.8),
  mutinfo_estimator = function(x, y) {
    -0.5 * log(1 - stats::cor(x, y, method = "spearman"))
  }
)
print("1")
set.seed(1)
sensitivity_intervals_sprho_subset = sensitivity_intervals_Dvine(
  fitted_model = best_fitted_model,
  sens_results = sens_results_sprho,
  mutinfo_estimator = function(x, y) {
    -0.5 * log(1 - stats::cor(x, y, method = "spearman"))
  },
  n_prec = n_prec,
  B = B,
  ncores = ncores
)
print(Sys.time() - a)

# Save Results ------------------------------------------------------------

# The results of the sensitivity analysis are saved to a file. These results are
# analyzed in a separate file.
saveRDS(
  object = sens_results_sicc,
  file = "results/sensitivity-analysis-results-main-sicc.rds"
)
saveRDS(
  object = sensitivity_intervals_Rh_subset,
  file ="results/sensitivity-intervals-Rh-subset.rds"
)
saveRDS(
  object = sensitivity_intervals_sprho_full,
  file ="results/sensitivity-intervals-sprho-full.rds"
)

saveRDS(
  object = sens_results_sprho,
  file = "results/sensitivity-analysis-results-main-sprho.rds"
)
saveRDS(
  object = sensitivity_intervals_sprho_subset,
  file ="results/sensitivity-intervals-sprho-subset.rds"
)
