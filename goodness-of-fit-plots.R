# Preliminaries ---------------------------------------------------------
library(Surrogate)
library(tidyverse)
library(survminer)
library(survival)
library(latex2exp)
# Size parameter for saving plots to disk.
single_width = 61
double_width = 14 / 2.54
between_width = single_width * 1.5 / 25.4
single_height = 55
double_height = 12.8 / 2.54
between_height = single_height * 1.5 / 25.4
res = 600

save_to_main = "paper-figures-tables/main-text/"
save_to_appendix = "paper-figures-tables/appendix/goodness-of-fit/"

# Load fitted models.
fitted_models = readRDS("results/fitted-models.rds")
best_fitted_model = readRDS("results/best-fitted-model.rds")

# Evaluate goodness of fit of the selected model.
grid = seq(from = 1, to = 250, length.out = 400)
plot(best_fitted_model, grid = grid)
# The plot() method automatically runs all GoF plots. In order to save them
# one-by-one, we call the functions that produce only single plots one-by-one.
# However, the plots for the main text are produced manually using ggplot.

## Plots for Main text ----------------------------------------------------

# Function to make the GoF plots using ggplot2. This produces more beautiful
# plots than those produced with plot(). This function also automatically saves
# the plots to the correct folder.
ggplot_gof_plots = function(fitted_submodel, knots, knott, treat) {
  if (treat == 0) {
    fitted_submodel = best_fitted_model$fit_0
    knots = best_fitted_model$knots0
    knott = best_fitted_model$knott0
  }
  else if (treat == 1) {
    fitted_submodel = best_fitted_model$fit_1
    knots = best_fitted_model$knots1
    knott = best_fitted_model$knott1
  }

  para = fitted_submodel$estimate
  surv_joint = function(x) {
    exp(
      Surrogate:::survival_survival_loglik(
        para = para,
        X = x,
        delta_X = 0,
        Y = x,
        delta_Y = 0,
        copula_family = best_fitted_model$copula_family[treat + 1],
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
        treat
    ),
    ylab = TeX(paste0("$P(S_{", treat, "} > t)$")),
    xlab = "t (weeks)",
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
    filename = paste0(save_to_main, "marginal-gof-s", treat, ".pdf"),
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
        treat
    ),
    ylab = TeX(paste0("$P(T_{", treat, "} > t)$")),
    xlab = "t (weeks)",
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
    filename = paste0(save_to_main, "marginal-gof-t", treat, ".pdf"),
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
    treated = treat
  )
  fit_gam = mgcv::gam(
    y ~ s(x),
    family = stats::quasi(link = "log",
                          variance = "mu"),
    data = data.frame(x = selected_data$Surv[selected_data$Treat ==
                                               treat], y = selected_data$Pfs[selected_data$Treat ==
                                                                           treat])
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
                                                          treat], y = selected_data$Pfs[selected_data$Treat ==
                                                                                      treat]),
               alpha = 0.1) +
    geom_line(aes(x = grid, y = model_cond_means),
              color = "red",
              linetype = "dashed") +
    coord_cartesian(ylim = c(0, max(selected_data$Surv))) +
    xlab("t (weeks)") +
    ylab(TeX(paste0("$E(S_{", treat, "} | T_{", treat, "} = t, S_{", treat, "} < T_{", treat, "})$"))) +
    theme_bw()
  ggsave(
    filename = paste0(save_to_main, "mean-S-before-T-gof", treat,".pdf"),
    width = single_width,
    height = single_height,
    units = "mm",
    device = "pdf"
  )

  selected_data = best_fitted_model$data[best_fitted_model$data$SurvInd ==
                                           1,]
  model_prob = 1 - sapply(X = grid, FUN = Surrogate:::prob_progression_before_dying,
                          fitted_model = best_fitted_model, treated = treat)
  fit_gam = mgcv::gam(
    y ~ s(x),
    family = stats::binomial(),
    data = data.frame(x = selected_data$Surv[selected_data$Treat ==
                                               treat],
                      y = 1 - selected_data$PfsInd[selected_data$Treat ==
                                                     treat])
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
                                                          treat],
                                 y = 1 - selected_data$PfsInd[selected_data$Treat ==
                                                                treat]),
               alpha = 0.1) +
    geom_line(aes(x = grid, y = model_prob),
              color = "red",
              linetype = "dashed") +
    xlab("t (weeks)") +
    ylab(TeX(paste0("$P(S_{", treat, "} = T_{", treat, "} | T_{", treat, "} = t)$"))) +
    theme_bw()
  ggsave(
    filename = paste0(save_to_main, "prob-dying-gof", treat, ".pdf"),
    width = single_width,
    height = single_height,
    units = "mm",
    device = "pdf"
  )
}

ggplot_gof_plots(best_fitted_model$fit_0, best_fitted_model$knots0, best_fitted_model$knott0, 0)
ggplot_gof_plots(best_fitted_model$fit_1, best_fitted_model$knots1, best_fitted_model$knott1, 1)

## Plots for Appendix -----------------------------------------------------

# Function to produce and save GoF plots for the Web Appendix.
gof_plots_appendix = function(copula_family, fitted_model) {
  pdf(
    file = paste0(save_to_appendix, copula_family, "/", "marginal-gof-s0.pdf"),
    width = between_width,
    height = between_height
  )
  par(mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
  marginal_gof_scr_S_plot(
    fitted_model = fitted_model,
    grid = grid,
    treated = 0,
    # main = "Survival Function S",
    xlab = "t (weeks)",
    ylab = TeX("$P(S_{0} > t)$"),
    xlim = c(0, 250)
  )
  dev.off()
  pdf(file = paste0(save_to_appendix, copula_family, "/", "marginal-gof-s1.pdf"), width = between_width, height = between_height)
  par(mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
  marginal_gof_scr_S_plot(
    fitted_model = fitted_model,
    grid = grid,
    treated = 1,
    # main = "Survival Function S - Active Treatment",
    xlab = "t (weeks)",
    ylab = TeX("$P(S_{1} > t)$"),
    xlim = c(0, 250)
  )
  dev.off()
  pdf(file = paste0(save_to_appendix, copula_family, "/", "marginal-gof-t0.pdf"), width = between_width, height = between_height)
  par(mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
  marginal_gof_scr_T_plot(
    fitted_model = fitted_model,
    grid = grid,
    treated = 0,
    # main = "Survival Function T - Control Treatment",
    xlab = "t (weeks)",
    ylab = TeX("$P(T_{0} > t)$"),
    xlim = c(0, 250)
  )
  dev.off()
  pdf(file = paste0(save_to_appendix, copula_family, "/", "marginal-gof-t1.pdf"), width = between_width, height = between_height)
  par(mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
  marginal_gof_scr_T_plot(
    fitted_model = fitted_model,
    grid = grid,
    treated = 1,
    # main = "Survival Function T - Active Treatment",
    xlab = "t (weeks)",
    ylab = TeX("$P(T_{1} > t)$"),
    xlim = c(0, 250)
  )
  dev.off()

  pdf(file = paste0(save_to_appendix, copula_family, "/", "mean-S-before-T-gof0.pdf"), width = between_width, height = between_height)
  par(mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
  mean_S_before_T_plot_scr(
    fitted_model = fitted_model,
    grid = grid,
    treated = 0,
    xlab = "t (weeks)",
    ylab = TeX("$E(S_{0} | T_{0} = t, S_{0} < T_{0})"),
    col = "gray",
    # main = "Control Treatment",
    xlim = c(0, 250)
  )
  dev.off()
  pdf(file = paste0(save_to_appendix, copula_family, "/", "mean-S-before-T-gof1.pdf"), width = between_width, height = between_height)
  par(mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
  mean_S_before_T_plot_scr(
    fitted_model = fitted_model,
    grid = grid,
    treated = 1,
    xlab = "t (weeks)",
    ylab = TeX("$E(S_{1} | T_{1} = t, S_{1} < T_{1})"),
    col = "gray",
    # main = "Active Treatment",
    xlim = c(0, 250)
  )
  dev.off()
  pdf(file = paste0(save_to_appendix, copula_family, "/", "prob-dying-gof0.pdf"), width = between_width, height = between_height)
  par(mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
  prob_dying_without_progression_plot(
    fitted_model = fitted_model,
    grid = grid,
    treated = 0,
    xlab = "t (weeks)",
    ylab = TeX("$P(S_{0} = T_{0} | T_{0} = t)$"),
    col = "gray",
    # main = "Control Treatment",
    xlim = c(0, 250)
  )
  dev.off()
  pdf(file = paste0(save_to_appendix, copula_family, "/", "prob-dying-gof1.pdf"), width = between_width, height = between_height)
  par(mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
  prob_dying_without_progression_plot(
    fitted_model = fitted_model,
    grid = grid,
    treated = 1,
    xlab = "t (weeks)",
    ylab = TeX("$P(S_{1} = T_{1} | T_{1} = t)$"),
    col = "gray",
    # main = "Active Treatment",
    xlim = c(0, 250)
  )
  dev.off()
}

# Select best model for each copula except the Gaussian copula (for which the plots
# are produced with ggplot 2).
fitted_models %>%
  filter(copula != "gaussian") %>%
  group_by(copula) %>%
  slice_min(AIC) %>%
  summarize(gof_plots_appendix(copula, fitted_model[[1]]))


