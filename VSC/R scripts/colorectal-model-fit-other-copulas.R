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
# model is fitted. The Frank, Gumbel, and Clayton models with two internal knots
# are fitted next and the sensitivity analysis performed.
fitted_model = fit_model_SurvSurv(data = data,
                                  copula_family = "gaussian",
                                  nknots = 2)

set.seed(1)
sens_results_frank = sensitivity_analysis_SurvSurv_copula(
  fitted_model = fitted_model,
  copula_family2 = "frank",
  degrees = 0,
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

set.seed(1)
sens_results_gumbel = sensitivity_analysis_SurvSurv_copula(
  fitted_model = fitted_model,
  copula_family2 = "gumbel",
  degrees = 0,
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

set.seed(1)
sens_results_clayton = sensitivity_analysis_SurvSurv_copula(
  fitted_model = fitted_model,
  copula_family2 = "clayton",
  degrees = 0,
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

# The results of the sensitivity analyses with different copulas are combined
# into a single data frame. The combined results of the sensitivity analysis are
# saved to a file. These results are analyzed in a separate file.
sens_results_combined =
  dplyr::bind_rows(
    sens_results_frank %>% mutate(copula = "Frank"),
    sens_results_gumbel %>% mutate(copula = "Gumbel"),
    sens_results_clayton %>% mutate(copula = "Clayton")
  )
readr::write_csv(
  x = sens_results_combined,
  file = paste0(saveto, "/sensitivity-analysis-results-other-copulas.csv")
)




