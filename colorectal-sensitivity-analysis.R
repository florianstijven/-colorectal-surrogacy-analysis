# Setup -------------------------------------------------------------------
# Load the required packages
library(Surrogate)
library(tidyverse)
# Load the best fitted model.
best_fitted_model = readRDS("results/best-fitted-model.rds")


# Sensitivity Analysis ----------------------------------------------------

# The sensitivity analysis is implemented in this file, but the results of the
# sensitivity analysis are saved and processed elsewhere. This is done because
# the sensitivity analysis is computer intensive and not interactive; whereas
# processing the results is not computer intensive, but interactive.

# Run the sensitivity analysis with the main assumptions.
set.seed(1)
a = Sys.time()
sens_results = sensitivity_analysis_SurvSurv_copula(
  fitted_model = best_fitted_model,
  n_sim = 5e3,
  n_prec = 1e4,
  ncores = 10,
  marg_association = TRUE,
  cond_ind = TRUE,
  composite = TRUE,
  degrees = 0,
  copula_family2 = "gaussian",
  lower = c(0.5, 0, 0, 0.15),
  upper = c(0.95, 0, 0, 0.8)
)
print(Sys.time() - a)


# Save Results ------------------------------------------------------------

# The results of the sensitivity analysis are saved to a file. These results are
# analyzed in a separate file.
saveRDS(
  object = sens_results,
  file = "results/sensitivity-analysis-results-main.rds"
)
