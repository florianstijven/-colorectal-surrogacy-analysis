#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)

# Setup -------------------------------------------------------------------

Sys.setenv(TZ='Europe/Brussels')
# Sandboxing takes a significant amount of time on the HPC and is therefore
# disabled.
options(renv.config.sandbox.enabled = FALSE)
ncores = as.integer(args[1])

library(Surrogate)
library(dplyr)
library(tidyr)
library(copula)
library(FNN)

# We need the best fitted vine copula model.
best_fitted_model = readRDS("results/best-fitted-model.rds")

# We define all different scenarios for the set of sensitivity analysis in the
# scenarios_tbl.
sensitivity_ranges = tibble(
  ranges = list(list(
    lower = c(0.5, 0, 0, 0.15),
    upper = c(0.95, 0, 0, 0.8)
  ), list(
    lower = c(0.25, 0, 0, 0.05),
    upper = c(0.95, 0.3, 0.3, 0.8)
  )),
  range_class = c("Main Assumptions", "Relaxed Assumptions"),
  cond_ind = c(TRUE, FALSE)
)
# We consider all combinations of parameter ranges, unidentifiable copula
# families, and ICA (Spearman's rho or SICC).
scenarios_tbl = expand_grid(
  sensitivity_ranges,
  copula_family = c("gaussian", "frank", "gumbel", "clayton"),
  ICA_type = c("SICC", "Spearman's correlation")
)
# The SICC can be replaced with any measure by replacing the mutual information
# estimator with an estimator of -0.5 * log(1 - measure).
scenarios_tbl = scenarios_tbl %>%
  mutate(mutinfo_estimator = list(function(x, y) {
    -0.5 * log(1 - stats::cor(x, y, method = "spearman"))
  }))
scenarios_tbl$mutinfo_estimator[scenarios_tbl$ICA_type == "SICC"] = list(NULL)

# Sensitivity Analysis ----------------------------------------------------

# We use a wrapper function for the sensitivity analysis such that we set the
# same seed for each different version of the sensitivity analysis.
wrapper_sensitivity_analysis = function(cond_ind, copula_family, lower, upper, mutinfo_estimator) {
  set.seed(1)
  sensitivity_analysis_SurvSurv_copula(
    fitted_model = best_fitted_model,
    n_sim = 100,
    n_prec = 5000,
    ncores = ncores,
    marg_association = TRUE,
    cond_ind = cond_ind,
    composite = TRUE,
    copula_family2 = copula_family,
    degrees = 0,
    lower = ranges$lower,
    upper = ranges$upper,
    mutinfo_estimator = mutinfo_estimator
  )
}

# The sensitivity analysis is implemented in this file, but the results of the
# sensitivity analysis are saved into an .RData file and processed elsewhere.
# This is done because the sensitivity analysis is computer intensive and not
# interactive; whereas processing the results is not computer intensive, but
# interactive.
a = Sys.time()
sens_results_tbl = scenarios_tbl %>%
  rowwise(everything()) %>%
  summarize(
    sens_results = list(sensitivity_analysis_SurvSurv_copula(
      fitted_model = best_fitted_model,
      n_sim = 100,
      n_prec = 5000,
      ncores = ncores,
      marg_association = TRUE,
      cond_ind = cond_ind,
      composite = TRUE,
      copula_family2 = copula_family,
      degrees = 0,
      lower = ranges$lower,
      upper = ranges$upper,
      mutinfo_estimator = mutinfo_estimator
    ))
  )
print(Sys.time() - a)


# Saving Results ----------------------------------------------------------

# The results of the sensitivity analysis are saved to a file. These results are
# analyzed in a separate file.
saveRDS(
  object = sens_results_tbl,
  file = "results/sensitivity-analysis-results-relaxed.rds"
)
