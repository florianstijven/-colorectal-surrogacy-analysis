# Setup -------------------------------------------------------------------

# Load the required packages
library(Surrogate)
library(tidyverse)
# Size parameter for saving plots to disk.
single_width = 9 / 2.54
double_width = 14 / 2.54
between_width = single_width * 1.25
single_height = 8.2 / 2.54
double_height = 12.8 / 2.54
between_height = single_height * 1.25
res = 600

save_to_main = "paper-figures-tables/main-text/"
save_to_appendix = "paper-figures-tables/appendix/"

# Put the data in correct format: (x_i, y_i, z_i, \delta_i^X, \delta_i^Y). This
# is also explained in the help files of Surrogate::fit_model_SurvSurv.
data = read.csv("COLO_DAT.csv")
data = data %>%
  select(PFSTIME, SURVTIME, TREAT, PFS_IND, SURV_IND)
# In the original data, the surrogate variable is the composite progression-free
# survival (PFS). As explained in the paper, we will model time-to progression
# (TTP) instead, where TTP is dependently censored by overall survival (OS).
# Therefore, PFS is recoded into TTP in these data.
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


# Model Fitting -----------------------------------------------------------

# After ensuring that the data are in the correct format, the survival-survival
# model is fitted. Models are fitted for different combinations of the number of
# internal knots and the parametric copula families.

# Four parametric copula families are considered.
possible_copulas = c("gaussian", "clayton", "frank", "gumbel")
# 2 through 5 internal knots are considered. The same number of internal knots
# are considered for each potential outcome, although this can be relaxed.
possible_nknots = 2:5
# Construct a tibble with all possible combinations of the number of internal
# knots and the parametric copula families.
model_combinations = expand_grid(copula = possible_copulas,
                                 nknots = possible_nknots)
# Models for all the above combinations are fitted. The fitted models are saved
# in a column of the tibble.
fitted_models = model_combinations %>%
  mutate(fitted_model = purrr::map2(
    .x = copula,
    .y = nknots,
    .f = function(.x, .y) {
      fit_model_SurvSurv(data, .x, .y)
    }
  ))

# For all fitted models, goodness-of-fit measures are computed and saved into
# the same tibble. The maximized loglikelihood of the entire identifiable model
# is the sum of the loglikelihoods of the two fitted submodels.
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
    AIC = -2 * LogLik + 2 * df
  ) %>%
  arrange(AIC)

# Print summary of all fitted models order from lowest to largest AIC. A lower
# AIC corresponds to a better fit.
sink(file = paste0(save_to_appendix, "fitted-models.txt")) # Open connection to .txt file to print output to
cat("Table of fitted models:\n\n")
print(fitted_models %>%
        mutate(
          LogLik = num(LogLik, digits = 2),
          AIC = num(AIC, digits = 2)
        ))
sink()


# The best fitting model, in terms of AIC, is the Gaussian copula model with 2
# internal knots. This model is extracted from the list of fitted models and
# named best_fitted_model. Since the fitted_models tibble is already sorted on
# AIC, the best model is the first one.
sink(file = "results/best-fitted-model-summary.txt")
best_fitted_model = fitted_models$fitted_model[[1]]
# Print summary of the selected model.
cat("\nBest fitted model:\n\n")
best_fitted_model
sink() # Close connection to .txt file.


# Goodness of Fit ---------------------------------------------------------

# Evaluate goodness of fit of the selected model.
grid = seq(from = 1, to = 250, length.out = 400)
plot(best_fitted_model, grid = grid)
# The plot() method automatically runs all GoF plots. In order to save them
# one-by-one. We call the functions that produce only single plots one-by-one.

pdf(file = paste0(save_to_main, "marginal-gof-s0.pdf"), width = between_width, height = between_height)
marginal_gof_scr_S_plot(
  fitted_model = best_fitted_model,
  grid = grid,
  treated = 0,
  main = "Survival Function S - Control Treatment",
  xlab = "t (months)",
  ylab = "P(S > t)",
  xlim = c(0, 250)
)
dev.off()
pdf(file = paste0(save_to_appendix, "marginal-gof-s1.pdf"), width = between_width, height = between_height)
marginal_gof_scr_S_plot(
  fitted_model = best_fitted_model,
  grid = grid,
  treated = 1,
  main = "Survival Function S - Active Treatment",
  xlab = "t (months)",
  ylab = "P(S > t)",
  xlim = c(0, 250)
)
dev.off()
pdf(file = paste0(save_to_main, "marginal-gof-t0.pdf"), width = between_width, height = between_height)
marginal_gof_scr_T_plot(
  fitted_model = best_fitted_model,
  grid = grid,
  treated = 0,
  main = "Survival Function T - Control Treatment",
  xlab = "t (months)",
  ylab = "P(T > t)",
  xlim = c(0, 250)
)
dev.off()
pdf(file = paste0(save_to_appendix, "marginal-gof-t1.pdf"), width = between_width, height = between_height)
marginal_gof_scr_T_plot(
  fitted_model = best_fitted_model,
  grid = grid,
  treated = 1,
  main = "Survival Function T - Active Treatment",
  xlab = "t (months)",
  ylab = "P(T > t)",
  xlim = c(0, 250)
)
dev.off()

pdf(file = paste0(save_to_main, "mean-S-before-T-gof0.pdf"), width = between_width, height = between_height)
mean_S_before_T_plot_scr(
  fitted_model = best_fitted_model,
  grid = grid,
  treated = 0,
  xlab = "t (months)",
  ylab = "E(S | T = t, S < T)",
  col = "gray",
  main = "Control Treatment",
  xlim = c(0, 250)
)
dev.off()
pdf(file = paste0(save_to_appendix, "mean-S-before-T-gof1.pdf"), width = between_width, height = between_height)
mean_S_before_T_plot_scr(
  fitted_model = best_fitted_model,
  grid = grid,
  treated = 1,
  xlab = "t (months)",
  ylab = "E(S | T = t, S < T)",
  col = "gray",
  main = "Active Treatment",
  xlim = c(0, 250)
)
dev.off()
pdf(file = paste0(save_to_main, "prob-dying-gof0.pdf"), width = between_width, height = between_height)
prob_dying_without_progression_plot(
  fitted_model = best_fitted_model,
  grid = grid,
  treated = 0,
  xlab = "t (months)",
  ylab = "P(S = T | T = t)",
  col = "gray",
  main = "Control Treatment",
  xlim = c(0, 250)
)
dev.off()
pdf(file = paste0(save_to_appendix, "prob-dying-gof1.pdf"), width = between_width, height = between_height)
prob_dying_without_progression_plot(
  fitted_model = best_fitted_model,
  grid = grid,
  treated = 1,
  xlab = "t (months)",
  ylab = "P(S = T | T = t)",
  col = "gray",
  main = "Active Treatment",
  xlim = c(0, 250)
)
dev.off()


# Saving Results ----------------------------------------------------------

saveRDS(fitted_models, file = "results/fitted-models.rds")
saveRDS(best_fitted_model, file = "results/best-fitted-model.rds")
