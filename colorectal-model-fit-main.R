# Setup -------------------------------------------------------------------

# Load the required packages
library(Surrogate)
library(tidyverse)
library(survminer)
library(survival)
# Size parameter for saving plots to disk.
single_width = 61
double_width = 14 / 2.54
between_width = single_width * 1.25
single_height = 55
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

sink(file = paste0(save_to_main, "data-descriptive-statistics.txt"))
# Compute the number (and proportion) of censoring for overall survival.
data %>%
  summarize(
    "Total number of observations" = n(),
    "Number of censored obervations for OS" = sum(SURV_IND == 0),
    "Number of OS events" = sum(SURV_IND),
    "Number of censored observations for PFS" = sum(PFS_IND == 0),
    "Number of PFS events" = sum(PFS_IND),
    "Proportion of censored observations for OS" = mean(SURV_IND == 0)
  )
cat("\n")
# Estimate the median PFS and OS durations.
cat("Median PFS time in weeks (95% CI)\n")
surv_median(survfit(Surv(PFSTIME, PFS_IND) ~ 1, data))
cat("\nMedian OS time in weeks (95% CI)\n")
surv_median(survfit(Surv(SURVTIME, SURV_IND) ~ 1, data))
sink()


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
# one-by-one, we call the functions that produce only single plots one-by-one.
# However, the plots for the main text are produced manually using ggplot.

## Plots for Main text ----------------------------------------------------
fitted_submodel = best_fitted_model$fit_0
knots = best_fitted_model$knots0
knott = best_fitted_model$knott0

para = fitted_submodel$estimate
surv_joint = function(x) {
  exp(
    Surrogate:::survival_survival_loglik(
      para = para,
      X = x,
      delta_X = 0,
      Y = x,
      delta_Y = 0,
      copula_family = best_fitted_model$copula_family[1],
      knotsx = knots,
      knotsy = knott
    )
  )
}
Surv_probs = sapply(X = grid, FUN = surv_joint)
ggsurvplot(
  survival::survfit(
    survival::Surv(
      Pfs,
      Ind
    ) ~
      1,
    data = best_fitted_model$data %>%
      mutate(Ind = pmax(PfsInd, SurvInd)),
    subset = best_fitted_model$data$Treat ==
      0
  ),
  ylab = "P(S > t)",
  xlab = "t (months)",
  xlim = c(0, 250),
  color = "black",
  censor = FALSE
) %>% `[[`(1) +
  geom_line(aes(x = grid, y = Surv_probs),
            data = data.frame(grid = grid, Surv_probs), color = "red", linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250)) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(
  filename = paste0(save_to_main, "marginal-gof-s0.pdf"),
  width = single_width,
  height = single_height,
  units = "mm",
  device = "pdf"
)

ks = length(knots)
kt = length(knott)
para_t = para[(1 + ks):(ks + kt)]
surv_t = function(x) {
  flexsurv::psurvspline(q = x, gamma = para_t, knots = knott)
}
Surv_probs = 1 - surv_t(grid)
ggsurvplot(
  survival::survfit(
    survival::Surv(
      Surv,
      SurvInd
    ) ~
      1,
    data = best_fitted_model$data,
    subset = best_fitted_model$data$Treat ==
      0
  ),
  ylab = "P(T > t)",
  xlab = "t (months)",
  xlim = c(0, 250),
  color = "black",
  censor = FALSE
) %>% `[[`(1) +
  geom_line(aes(x = grid, y = Surv_probs),
            data = data.frame(grid = grid, Surv_probs), color = "red", linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250)) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(
  filename = paste0(save_to_main, "marginal-gof-t0.pdf"),
  width = single_width,
  height = single_height,
  units = "mm",
  device = "pdf"
)


selected_data = best_fitted_model$data[best_fitted_model$data$PfsInd ==
                                         1 &
                                         best_fitted_model$data$SurvInd == 1,]
model_cond_means = sapply(
  X = grid,
  FUN = Surrogate:::mean_S_before_T,
  fitted_model = best_fitted_model,
  treated = 0
)
fit_gam = mgcv::gam(
  y ~ s(x),
  family = stats::quasi(link = "log",
                        variance = "mu"),
  data = data.frame(x = selected_data$Surv[selected_data$Treat ==
                                             0], y = selected_data$Pfs[selected_data$Treat ==
                                                                               0])
)
predictions = mgcv::predict.gam(fit_gam,
                                newdata = data.frame(x = grid),
                                type = "link",
                                se.fit = TRUE)
data.frame(
  grid,
  pred = exp(predictions$fit),
  upper_ci = exp(predictions$fit + 1.96 *
                   predictions$se.fit),
  lower_ci = exp(predictions$fit - 1.96 *
                   predictions$se.fit),
  model_cond_means
) %>%
  ggplot() +
  geom_line(aes(x = grid, y = pred)) +
  geom_ribbon(aes(x = grid, ymin = lower_ci, ymax = upper_ci), alpha = 0.3) +
  geom_point(aes(x = x, y = y),
             data = data.frame(x = selected_data$Surv[selected_data$Treat ==
                                                        0], y = selected_data$Pfs[selected_data$Treat ==
                                                                                    0]),
             alpha = 0.1) +
  geom_line(aes(x = grid, y = model_cond_means),
            color = "red",
            linetype = "dashed") +
  coord_cartesian(ylim = c(0, max(selected_data$Surv))) +
  xlab("t (months)") +
  ylab("E(S | T = t, S < T)") +
  theme_bw()
ggsave(
  filename = paste0(save_to_main, "mean-S-before-T-gof0.pdf"),
  width = single_width,
  height = single_height,
  units = "mm",
  device = "pdf"
)

selected_data = best_fitted_model$data[best_fitted_model$data$SurvInd ==
                                    1,]
model_prob = 1 - sapply(X = grid, FUN = Surrogate:::prob_progression_before_dying,
                        fitted_model = best_fitted_model, treated = 0)
fit_gam = mgcv::gam(
  y ~ s(x),
  family = stats::binomial(),
  data = data.frame(x = selected_data$Surv[selected_data$Treat ==
                                             0],
                    y = 1 - selected_data$PfsInd[selected_data$Treat ==
                                                   0])
)
expit = function(x)
  1 / (1 + exp(-x))
predictions = mgcv::predict.gam(fit_gam,
                                newdata = data.frame(x = grid),
                                type = "link",
                                se.fit = TRUE)
data.frame(
  grid,
  pred = expit(predictions$fit),
  upper_ci = expit(predictions$fit + 1.96 *
                   predictions$se.fit),
  lower_ci = expit(predictions$fit - 1.96 *
                   predictions$se.fit),
  model_prob
) %>%
  ggplot() +
  geom_line(aes(x = grid, y = pred)) +
  geom_ribbon(aes(x = grid, ymin = lower_ci, ymax = upper_ci), alpha = 0.3) +
  geom_point(aes(x = x, y = y),
             data = data.frame(x = selected_data$Surv[selected_data$Treat ==
                                                        0],
                               y = 1 - selected_data$PfsInd[selected_data$Treat ==
                                                              0]),
             alpha = 0.1) +
  geom_line(aes(x = grid, y = model_prob),
            color = "red",
            linetype = "dashed") +
  xlab("t (months)") +
  ylab("P(S = T | T = t)") +
  theme_bw()
ggsave(
  filename = paste0(save_to_main, "prob-dying-gof0.pdf"),
  width = single_width,
  height = single_height,
  units = "mm",
  device = "pdf"
)


## Plots for Appendix -----------------------------------------------------

pdf(file = paste0(save_to_appendix, "marginal-gof-s0.pdf"), width = between_width, height = between_height)
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
pdf(file = paste0(save_to_appendix, "marginal-gof-t0.pdf"), width = between_width, height = between_height)
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

pdf(file = paste0(save_to_appendix, "mean-S-before-T-gof0.pdf"), width = between_width, height = between_height)
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
pdf(file = paste0(save_to_appendix, "prob-dying-gof0.pdf"), width = between_width, height = between_height)
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
