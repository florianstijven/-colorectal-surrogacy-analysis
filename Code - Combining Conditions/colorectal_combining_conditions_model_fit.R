library(Package)
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

# The best fitting model was determined in colorectal_model_fit.R Here, that
# model is again fitted because the fitted model object is needed for the
# sensitivity analysis function.
best_fitted_model = fit_model_SurvSurv(data = data,
                                       copula_family = "gaussian",
                                       nknots = 2)
# SENSITIVITY ANALYSIS

# The sensitivity analysis is implemented in this file, but the results of the
# sensitivity analysis are saved into an .RData file and processed elsewhere.
# This is done because the sensitivity analysis is computer intensive and not
# interactive; whereas processing the results is not computer intensive, but
# interactive.

# Sensitivity analysis without equality assumptions with 50k replications, but
# R^2_H is not computed to save some computer time. The essence of the results
# are transportable between the ICA as Spearman's rho and R^2_H.
sens_results_no_cond_ind = ica_SurvSurv_sens(
  fitted_model = best_fitted_model,
  n_sim = 5e4,
  n_prec = 5000,
  minfo_prec = 0,
  restr = TRUE,
  copula_family2 = best_fitted_model$copula_family,
  ncores = 10,
  get_marg_tau = TRUE,
  cond_ind = FALSE
)


# The results of the sensitivity analysis are saved to a file. These results are
# analyzed in a separate file.
save(sens_results_no_cond_ind,
     file = "Code - Combining Conditions/sensitivity_analysis_results_100k.RData")
