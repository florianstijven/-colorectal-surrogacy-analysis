#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Path to R libraries
.libPaths(new = "/data/leuven/351/vsc35197/R")
print(args)
Sys.setenv(TZ='Europe/Brussels')
ncores = as.integer(args[1])
saveto = args[2]
library(Surrogate)
library(dplyr)
library(purrr)
library(copula)
library(flexsurv)


# put data in correct format: (x_i, y_i, z_i, \delta_i^X, \delta_i^Y)
data = read.csv(
  paste0(saveto, "/COLO_DAT.csv")
)
data = data %>%
  select(PFSTIME, SURVTIME, TREAT, PFS_IND, SURV_IND)
# In the original data, the surrogate variable is the composite progression-free
# survival (PFS). As explained in the paper, we will model time-to progression
# (TTP) instead, where TTP is dependently censored by overall survival (OS).
# Therefore, PFS is recoded to TTP in these data.
data = data %>%
  mutate(PROGRESS_IND = ifelse((PFSTIME == SURVTIME) &
                                 (PFS_IND == 1) & (SURV_IND == 1),
                               0,
                               PFS_IND),
         PROGRESSTIME = PFSTIME)
data = data %>%
  select(PROGRESSTIME, SURVTIME, TREAT, PROGRESS_IND, SURV_IND)

# Compute the proportion of dependent censoring in both treatment groups.
data %>%
  group_by(TREAT) %>%
  summarize(prop_dependent_censoring =
              mean((PROGRESS_IND == 0) & (SURV_IND == 1) &
                     (PROGRESSTIME == SURVTIME)))

# After ensuring that the data are in the correct format, the survival-survival
# model is fitted. Models ar fitted for different combinations of the number of
# knots and the parametric copula families.

# Four parametric copula families are considered.
possible_copulas = c("gaussian", "clayton", "frank", "gumbel")
# 2 through 5 internal knots are considered.
possible_nknots = 2:5
# Construct a list with all possible combinations of the number of internal
# knots and the parametric copula families.
model_combinations = tidyr::expand_grid(copula = possible_copulas, nknots = possible_nknots)
# Models for all the above combinations are fitted and put into a list.
fitted_models = model_combinations %>%
  mutate(fitted_model = purrr::map2(
    .x = copula,
    .y = nknots,
    .f = function(.x, .y) {
      fit_model_SurvSurv(data, .x, .y)
    }
  ))

# For each fitted model, model fit measures are obtained. Those are summarized
# in a table which is sorted from lowest to highest AIC.
fitted_models = fitted_models %>%
  mutate(
    LogLik = purrr::map_dbl(
      .x = fitted_model,
      .f = function(.x) {
        logLik(.x$fit_0) +
          logLik(.x$fit_1)
      }
    ),
    df = 2 * (2 * (nknots + 2) + 1),
    AIC = (-2 * LogLik) + (2 * df)
  ) %>%
  arrange(AIC)

# Print summary of all fitted models order from lowest to largest AIC.
fitted_models

# The best fitting model, in terms of AIC, is the Gaussian copula model with 2
# internal knots. This model is extracted from the list of fitted models and
# named best_fitted_model.
best_fitted_model = fitted_models$fitted_model[[1]]
# Print summary of the selected model.
best_fitted_model

# The marginal goodness of fit is assessed by comparing the model-based with the
# KM-estimate of the survival curve of OS and PFS (not TTP!) in both treatment
# groups. This is implemented in the marginal_gof_scr() function.
grid = seq(0.01, 200, 0.1)
marginal_gof_scr(fitted_model = best_fitted_model,
                 data = data, grid = grid, time_unit = "weeks")


# SENSITIVITY ANALYSIS

# The sensitivity analysis is implemented in this file, but the results of the
# sensitivity analysis are saved into an .RData file and processed elsewhere.
# This is done because the sensitivity analysis is computer intensive and not
# interactive; whereas processing the results is not computer intensive, but
# interactive.
set.seed(1)
a = Sys.time()
sens_results = sensitivity_analysis_SurvSurv_copula(
  fitted_model = best_fitted_model,
  n_sim = 5e3,
  n_prec = 5000,
  minfo_prec = 5000,
  ncores = ncores,
  marg_association = TRUE,
  cond_ind = TRUE,
  composite = TRUE,
  lower = c(0.5, 0, 0, 0.15),
  upper = c(0.95, 0, 0, 0.8)
)
print(Sys.time() - a)

# The results of the sensitivity analysis are saved to a file. These results are
# analyzed in a separate file.
readr::write_csv(
  x = sens_results,
  file = paste0(saveto, "/sensitivity-analysis-results-main.csv")
)


