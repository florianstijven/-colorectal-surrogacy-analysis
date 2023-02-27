library(Surrogate)
library(tidyverse)


# put data in correct format: (x_i, y_i, z_i, \delta_i^X, \delta_i^Y)
data = read.csv("COLO_DAT.csv")
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
possible_copulas = list("gaussian", "clayton", "frank", "gumbel")
# 2 through 5 internal knots are considered.
possible_nknots = 2:5
# Construct a list with all possible combinations of the number of internal
# knots and the parametric copula families.
model_combinations = cross2(.x = possible_copulas, .y = possible_nknots)
# Models for all the above combinations are fitted and put into a list.
fitted_model_combinations = map2(.x = map(model_combinations, 1),
                                 .y = map(model_combinations, 2),
                                 .f = ~fit_model_SurvSurv(data = data, 
                                                 copula_family = .x,
                                                 nknots = .y))

# For each fitted model, model fit measures are obtained. Those are summarized
# in a table which is sorted from lowest to highest AIC.
map_df(.x = fitted_model_combinations, 
       .f = model_fit_measures) %>%
  mutate(copula_family = map_chr(model_combinations, 1),
         nknots = map_dbl(model_combinations, 2)) %>%
  arrange(AIC)
# The best fitting model is the Gaussian copula model with 2 internal knots.
# This model is extracted from the list of fitted models and named
# best_fitted_model.
best_fitted_model = keep(.x = fitted_model_combinations,
                         .p = ~.$copula_family == "gaussian" && length(.$knots0) - 2 == 2)[[1]]


# The marginal goodness of fit is assessed by comparing the model-based with the
# KM-estimate of the survival curve of OS and PFS (not TTP!) in both treatment
# groups. This is implemented in the marginal_gof_scr() function.
grid = seq(0.01, 300, 0.1)
png(filename = "Figures/GOF_colo_scr.png", width = 650, height = 520)
marginal_gof_scr(fitted_model = best_fitted_model,
                 data = data, grid = grid, time_unit = "weeks")
dev.off()

# SENSITIVITY ANALYSIS

# The sensitivity analysis is implemented in this file, but the results of the
# sensitivity analysis are saved into an .RData file and processed elsewhere.
# This is done because the sensitivity analysis is computer intensive and not
# interactive; whereas processing the results is not computer intensive, but
# interactive.

# Sensitivity analysis without equality assumptions
sens_results_no_cond_ind = ica_SurvSurv_sens(
  fitted_model = best_fitted_model,
  n_sim = 5000,
  n_prec = 5000,
  minfo_prec = 5000,
  restr = TRUE,
  copula_family2 = best_fitted_model$copula_family,
  ncores = 10,
  get_marg_tau = TRUE,
  cond_ind = FALSE
)

# Sensitivity analysis with conditional independence assumptions
sens_results_cond_ind = ica_SurvSurv_sens(
  fitted_model = best_fitted_model,
  n_sim = 5000,
  n_prec = 5000,
  minfo_prec = 5000,
  restr = TRUE,
  copula_family2 = best_fitted_model$copula_family,
  ncores = 10,
  get_marg_tau = TRUE,
  cond_ind = TRUE
)


# The results of the sensitivity analysis are saved to a file. These results are 
# analyzed in a separate file.
save(sens_results_no_cond_ind, sens_results_cond_ind,
     file = "sensitivity_analysis_results.RData")