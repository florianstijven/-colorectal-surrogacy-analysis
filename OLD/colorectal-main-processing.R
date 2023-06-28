set.seed(1)
# Load the required packages
library(tidyverse)
library(latex2exp)
library(GGally)
#specify options for saving the plots to files
width = 13.2
height = 9.9
pointsize = 10
res = 500

# Load results of the sensitivity analysis
sens_results = readr::read_csv(file = "sensitivity-analysis-results-main.csv")

# ------------------------
# NO EQUALITY RESTRICTIONS
# ------------------------

# Path to save results
path = "Figures/main"

# Check whether stochastic monotonicity (SM) holds in every simulated scenario.
sens_results %>%
  filter(sp_rho_t0s1 > 0,
         sp_rho_t0t1 > 0,
         sp_rho_s0s1 > 0,
         sp_rho_s0t1 > 0) %>%
  summarise(n())
# All 5000 simulated scenarios satisfy SM.

# Histogram for R_H under different inequality assumptions.
png(
  filename = paste0(path, "colo_results_ICA.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  ggplot(aes(x = ICA)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$R_h^2$")) +
  geom_histogram(
    fill = "gray",
    color = "black",
    binwidth = 0.025,
    boundary = 1
  ) +
  theme_bw()
dev.off()


# Histogram for Spearman's rho under different inequality assumptions.
png(
  filename = paste0(path, "colo_results_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  ggplot(aes(x = sp_rho)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$\\rho_s$")) +
  geom_histogram(
    fill = "gray",
    color = "black",
    binwidth = 0.05,
    boundary = 1
  ) +
  theme_bw()
dev.off()

#TABLES
# Table for R^2_h
sens_results_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions) %>%
  summarise(Range = paste0("[", round(min(icc), 3), ", ",
                           round(max(icc), 3), "]"),
            perc = paste0("[", round(quantile(icc, 0.01), 3), ", ",
                          round(quantile(icc, 0.99), 3), "]"),
            median = median(icc)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(assumptions) %>%
  knitr::kable(format = "latex")

# Table for Spearman's rho
sens_results_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions) %>%
  summarise(Range = paste0("[", round(min(sp_rho), 3), ", ",
                           round(max(sp_rho), 3), "]"),
            perc = paste0("[", round(quantile(sp_rho, 0.01), 3), ", ",
                          round(quantile(sp_rho, 0.99), 3), "]"),
            median = median(sp_rho)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(assumptions) %>%
  knitr::kable(format = "latex")

# Table for kendall's tau
sens_results_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions) %>%
  summarise(Range = paste0("[", round(min(kendall), 3), ", ",
                           round(max(kendall), 3), "]"),
            perc = paste0("[", round(quantile(kendall, 0.01), 3), ", ",
                          round(quantile(kendall, 0.99), 3), "]"),
            median = median(kendall)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(assumptions) %>%
  knitr::kable(format = "latex")


# Exploration of the results of the sensitivity analysis.

# Compute the smallest observable association in terms of Spearman's rho.
min_ca =
  sens_results %>%
  summarise(mean_0 = mean(sp_s0t0), mean_1 = mean(sp_s1t1)) %>%
  summarise(min = min(mean_0, mean_1)) %>%
  as.numeric()
min_ca # 0.684

# The R^2_h is plotted against the average cross association. Each replication
# is represented by a dot. The color indicates whether the weaker
# cross-associations assumption is satisfied in that replication of the
# sensitivity analysis. A vertical line indicates the value for the smallest
# observable association.
png(
  filename = paste0(path, "colo_results_exploration_cross_association_icc.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  mutate(ca = 0.5 * (sp_s0t1 + sp_s1t0),
         wa = 0.5 * (sp_s0s1 + sp_t0t1)) %>%
  mutate(WCA = ifelse(pmin(abs(sp_s0t0), abs(sp_s1t1)) > pmax(abs(sp_s0t1), abs(sp_s1t0)),
                      "Yes", "No")) %>%
  ggplot(aes(x = ca, y = icc, color = WCA)) +
  scale_x_continuous(
    limits = c(-1, 1),
    name = TeX("$\\frac{\\rho_{S_0, T_1} + \\rho_{S_1, T_0}}{2}$")
  ) +
  scale_y_continuous(name = TeX("$\\R_h^2"), limits = c(0, 1)) +
  geom_vline(xintercept = min_ca) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

# The R^2_h is plotted against the average within association, i.e., the average
# association of (S_0, S_1) and (T_0, T_1). Each replication is represented by a
# dot. The color indicates whether the stochastic monotonicity assumption is
# satisfied in that replication of the sensitivity analysis. A vertical line
# indicates a zero average within association.
png(
  filename = paste0(path, "colo_results_exploration_within_endpoint_association_icc.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  mutate(ca = 0.5 * (sp_s0t1 + sp_s1t0),
         wa = 0.5 * (sp_s0s1 + sp_t0t1)) %>%
  mutate(SM = ifelse(pmin(sp_s0s1, sp_t0t1) > 0,
                               "Yes", "No")) %>%
  ggplot(aes(x = wa, y = icc, color = SM)) +
  scale_x_continuous(
    limits = c(-1, 1),
    name = TeX("$\\frac{\\rho_{S_0S_1} + \\rho_{T_0T_1}}{2}$")
  ) +
  scale_y_continuous(name = TeX("$\\R_h^2"), limits = c(0, 1)) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

# The R^2_h is plotted against the proportion of patients harmed. This is the
# proportion of patients that have the surrogate event before the true endpoint
# event under experimental treatment, but have the true endpoint event before
# the surrogate event under control treatment.
png(
  filename = paste0(path, "colo_results_exploration_prop_harmed_icc.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  ggplot(aes(x = prop_harmed, y = icc)) +
  geom_point(alpha = 0.3) +
  scale_y_continuous(
    limits = c(0, 1),
    name = TeX("$R^2_h$")
  ) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

# The R^2_h is plotted against the log odds ratio of the principal strata of
# diseased/non-diseased in both treatment groups. A patient is 'diseased' under
# treatment Z, if the surrogate event happens before the true event.
png(
  filename = paste0(path, "colo_results_exploration_logodds_icc.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  ggplot(aes(x = log((prop_always * prop_never) / (prop_harmed * prop_protected)
  ), y = icc)) +
  geom_point(alpha = 0.3) +
  scale_y_continuous(limits = c(0, 1),
                     name = TeX("$R^2_h$")) +
  scale_x_continuous(name = "log odds ratio of diseased status") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()


# The above explorative plots are repeated for the ICA quantified by Spearman's rho.
png(
  filename = paste0(path, "colo_results_exploration_cross_association_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  mutate(ca = 0.5 * (sp_s0t1 + sp_s1t0),
         wa = 0.5 * (sp_s0s1 + sp_t0t1)) %>%
  mutate(WCA = ifelse(pmin(abs(sp_s0t0), abs(sp_s1t1)) > pmax(abs(sp_s0t1), abs(sp_s1t0)),
                      "Yes", "No")) %>%
  ggplot(aes(x = ca, y = sp_rho, color = WCA)) +
  scale_x_continuous(
    limits = c(-1, 1),
    name = TeX("$\\frac{\\rho_{S_0, T_1} + \\rho_{S_1, T_0}}{2}$")
  ) +
  scale_y_continuous(name = TeX("$\\rho_S"), limits = c(-1, 1)) +
  geom_vline(xintercept = min_ca) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

png(
  filename = paste0(path, "colo_results_exploration_within_endpoint_association_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  mutate(ca = 0.5 * (sp_s0t1 + sp_s1t0),
         wa = 0.5 * (sp_s0s1 + sp_t0t1)) %>%
  mutate(SM = ifelse(pmin(sp_s0s1, sp_t0t1) > 0,
                               "Yes", "No")) %>%
  ggplot(aes(x = wa, y = sp_rho, color = SM)) +
  scale_x_continuous(
    limits = c(-1, 1),
    name = TeX("$\\frac{\\rho_{S_0S_1} + \\rho_{T_0T_1}}{2}$")
  ) +
  scale_y_continuous(name = TeX("$\\rho_s"), limits = c(-1, 1)) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

png(
  filename = paste0(path, "colo_results_exploration_prop_harmed_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  ggplot(aes(x = prop_harmed, y = sp_rho)) +
  geom_point(alpha = 0.3) +
  scale_y_continuous(limits = c(-1, 1),
                     name = TeX("$\\rho_s$")) +
  scale_x_continuous(name = "Proportion Harmed") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

png(
  filename = paste0(path, "colo_results_exploration_logodds_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  ggplot(aes(x = log((prop_always * prop_never) / (prop_harmed * prop_protected)
  ), y = sp_rho)) +
  geom_point(alpha = 0.3) +
  scale_y_continuous(limits = c(-1, 1),
                     name = TeX("$\\rho_s$")) +
  scale_x_continuous(name = "log odds ratio of diseased status") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()


# -----------------------------------------------
# EQUALITY RESTRICTIONS: CONDITIONAL INDEPENDENCE
# -----------------------------------------------

# Path to save results
path = "Figures/Conditional Independence/"

# The mutual information is computed in the sensitivity analysis. The
# corresponding informational coefficient of correlation is computed.
sens_results = sens_results_cond_ind %>% mutate(icc = 1 - exp(-2*minfo))

# Construct an artificial data set to ease further plots.
# Rows that correspond to the no inequality assumptions scenario.
sens_results_helper = sens_results %>%
  mutate(assumptions = "No Assumptions") %>%
  ungroup()
# Rows that correspond to the stochastic monotonicity scenario.
sens_results_helper = sens_results_helper %>%
  bind_rows(
    sens_results %>% filter(
      sp_s0s1 > 0,
      sp_s0t0 > 0,
      sp_s0t1 > 0,
      sp_s1t0 > 0,
      sp_s1t1 > 0,
      sp_t0t1 > 0
    ) %>%
      mutate(assumptions = "SM") %>%
      ungroup()
  )
# Rows that correspond to the weaker cross-associations scenario.
sens_results_helper = sens_results_helper %>%
  bind_rows(
    sens_results %>%
      filter(pmin(abs(sp_s0t0), abs(sp_s1t1)) > pmax(abs(sp_s0t1), abs(sp_s1t0))) %>%
      mutate(assumptions = "WCA") %>%
      ungroup()
  )
# Rows that correspond to both weaker cross-associations and stochastic
# monotonicity.
sens_results_helper = sens_results_helper %>%
  bind_rows(
    sens_results %>%
      filter(
        pmin(abs(sp_s0t0), abs(sp_s1t1)) > pmax(abs(sp_s0t1), abs(sp_s1t0)),
        sp_s0s1 > 0,
        sp_s0t0 > 0,
        sp_s0t1 > 0,
        sp_s1t0 > 0,
        sp_s1t1 > 0,
        sp_t0t1 > 0
      ) %>%
      ungroup() %>%
      mutate(assumptions = "SM + WCA")
  )
# The number of rows in each scenario is computed and added to the helper data
# set.
sens_results_helper = left_join(
  x = sens_results_helper,
  y = sens_results_helper %>%
    group_by(assumptions) %>%
    summarise(n = n()) %>%
    ungroup(),
  by = c("assumptions")
)


# Histogram for R_H under different inequality assumptions.
png(
  filename = paste0(path, "colo_results_icc.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results_helper %>%
  mutate(assumptions = paste0(assumptions, " (n = ", n, ")")) %>%
  ggplot(aes(x = icc)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$R_h^2$")) +
  geom_histogram(
    fill = "gray",
    color = "black",
    binwidth = 0.025,
    boundary = 1
  ) +
  theme_bw() +
  facet_wrap(facets = ~ assumptions)
dev.off()


# Histogram for Spearman's rho under different inequality assumptions.
png(
  filename = paste0(path, "colo_results_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results_helper %>%
  mutate(assumptions = paste0(assumptions, " (n = ", n, ")")) %>%
  ggplot(aes(x = sp_rho)) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_x_continuous(name = TeX("$\\rho_s$")) +
  geom_histogram(
    fill = "gray",
    color = "black",
    binwidth = 0.05,
    boundary = 1
  ) +
  theme_bw() +
  facet_wrap(facets = ~ assumptions)
dev.off()

#TABLES
# Table for R^2_h
sens_results_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions) %>%
  summarise(Range = paste0("[", round(min(icc), 3), ", ",
                           round(max(icc), 3), "]"),
            perc = paste0("[", round(quantile(icc, 0.01), 3), ", ",
                          round(quantile(icc, 0.99), 3), "]"),
            median = median(icc)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(assumptions) %>%
  knitr::kable(format = "latex")

# Table for Spearman's rho
sens_results_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions) %>%
  summarise(Range = paste0("[", round(min(sp_rho), 3), ", ",
                           round(max(sp_rho), 3), "]"),
            perc = paste0("[", round(quantile(sp_rho, 0.01), 3), ", ",
                          round(quantile(sp_rho, 0.99), 3), "]"),
            median = median(sp_rho)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(assumptions) %>%
  knitr::kable(format = "latex")

# Table for kendall's tau
sens_results_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions) %>%
  summarise(Range = paste0("[", round(min(kendall), 3), ", ",
                           round(max(kendall), 3), "]"),
            perc = paste0("[", round(quantile(kendall, 0.01), 3), ", ",
                          round(quantile(kendall, 0.99), 3), "]"),
            median = median(kendall)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(assumptions) %>%
  knitr::kable(format = "latex")


# Exploration of the results of the sensitivity analysis.

# Compute the smallest observable association in terms of Spearman's rho.
min_ca =
  sens_results %>%
  summarise(mean_0 = mean(sp_s0t0), mean_1 = mean(sp_s1t1)) %>%
  summarise(min = min(mean_0, mean_1)) %>%
  as.numeric()
min_ca # 0.684

# The R^2_h is plotted against the average cross association. Each replication
# is represented by a dot. The color indicates whether the weaker
# cross-associations assumption is satisfied in that replication of the
# sensitivity analysis. A vertical line indicates the value for the smallest
# observable association.
png(
  filename = paste0(path, "colo_results_exploration_cross_association_icc.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  mutate(ca = 0.5 * (sp_s0t1 + sp_s1t0),
         wa = 0.5 * (sp_s0s1 + sp_t0t1)) %>%
  mutate(WCA = ifelse(pmin(abs(sp_s0t0), abs(sp_s1t1)) > pmax(abs(sp_s0t1), abs(sp_s1t0)),
                      "Yes", "No")) %>%
  ggplot(aes(x = ca, y = icc, color = WCA)) +
  scale_x_continuous(
    limits = c(-1, 1),
    name = TeX("$\\frac{\\rho_{S_0, T_1} + \\rho_{S_1, T_0}}{2}$")
  ) +
  scale_y_continuous(name = TeX("$\\R_h^2"), limits = c(0, 1)) +
  geom_vline(xintercept = min_ca) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

# The R^2_h is plotted against the average within association, i.e., the average
# association of (S_0, S_1) and (T_0, T_1). Each replication is represented by a
# dot. The color indicates whether the stochastic monotonicity assumption is
# satisfied in that replication of the sensitivity analysis. A vertical line
# indicates a zero average within association.
png(
  filename = paste0(path, "colo_results_exploration_within_endpoint_association_icc.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  mutate(ca = 0.5 * (sp_s0t1 + sp_s1t0),
         wa = 0.5 * (sp_s0s1 + sp_t0t1)) %>%
  mutate(SM = ifelse(pmin(sp_s0s1, sp_t0t1) > 0,
                               "Yes", "No")) %>%
  ggplot(aes(x = wa, y = icc, color = SM)) +
  scale_x_continuous(
    limits = c(-1, 1),
    name = TeX("$\\frac{\\rho_{S_0S_1} + \\rho_{T_0T_1}}{2}$")
  ) +
  scale_y_continuous(name = TeX("$\\R_h^2"), limits = c(0, 1)) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

# The R^2_h is plotted against the proportion of patients harmed. This is the
# proportion of patients that have the surrogate event before the true endpoint
# event under experimental treatment, but have the true endpoint event before
# the surrogate event under control treatment.
png(
  filename = paste0(path, "colo_results_exploration_prop_harmed_icc.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  ggplot(aes(x = prop_harmed, y = icc)) +
  geom_point(alpha = 0.3) +
  scale_y_continuous(
    limits = c(0, 1),
    name = TeX("$R^2_h$")
  ) +
  scale_x_continuous(name = "Proportion Harmed") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

# The R^2_h is plotted against the log odds ratio of the principal strata of
# diseased/non-diseased in both treatment groups. A patient is 'diseased' under
# treatment Z, if the surrogate event happens before the true event.
png(
  filename = paste0(path, "colo_results_exploration_logodds_icc.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  ggplot(aes(x = log((prop_always * prop_never) / (prop_harmed * prop_protected)
  ), y = icc)) +
  geom_point(alpha = 0.3) +
  scale_y_continuous(limits = c(0, 1),
                     name = TeX("$R^2_h$")) +
  scale_x_continuous(name = "log odds ratio of diseased status") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()


# The above explorative plots are repeated for the ICA quantified by Spearman's rho.
png(
  filename = paste0(path, "colo_results_exploration_cross_association_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  mutate(ca = 0.5 * (sp_s0t1 + sp_s1t0),
         wa = 0.5 * (sp_s0s1 + sp_t0t1)) %>%
  mutate(WCA = ifelse(pmin(abs(sp_s0t0), abs(sp_s1t1)) > pmax(abs(sp_s0t1), abs(sp_s1t0)),
                      "Yes", "No")) %>%
  ggplot(aes(x = ca, y = sp_rho, color = WCA)) +
  scale_x_continuous(
    limits = c(-1, 1),
    name = TeX("$\\frac{\\rho_{S_0, T_1} + \\rho_{S_1, T_0}}{2}$")
  ) +
  scale_y_continuous(name = TeX("$\\rho_S"), limits = c(-1, 1)) +
  geom_vline(xintercept = min_ca) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

png(
  filename = paste0(path, "colo_results_exploration_within_endpoint_association_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  mutate(ca = 0.5 * (sp_s0t1 + sp_s1t0),
         wa = 0.5 * (sp_s0s1 + sp_t0t1)) %>%
  mutate(SM = ifelse(pmin(sp_s0s1, sp_t0t1) > 0,
                               "Yes", "No")) %>%
  ggplot(aes(x = wa, y = sp_rho, color = SM)) +
  scale_x_continuous(
    limits = c(-1, 1),
    name = TeX("$\\frac{\\rho_{S_0S_1} + \\rho_{T_0T_1}}{2}$")
  ) +
  scale_y_continuous(name = TeX("$\\rho_s"), limits = c(-1, 1)) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

png(
  filename = paste0(path, "colo_results_exploration_prop_harmed_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  ggplot(aes(x = prop_harmed, y = sp_rho)) +
  geom_point(alpha = 0.3) +
  scale_y_continuous(limits = c(-1, 1),
                     name = TeX("$\\rho_s$")) +
  scale_x_continuous(name = "Proportion Harmed") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()

png(
  filename = paste0(path, "colo_results_exploration_logodds_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results %>%
  ggplot(aes(x = log((prop_always * prop_never) / (prop_harmed * prop_protected)
  ), y = sp_rho)) +
  geom_point(alpha = 0.3) +
  scale_y_continuous(limits = c(-1, 1),
                     name = TeX("$\\rho_s$")) +
  scale_x_continuous(name = "log odds ratio of diseased status") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
dev.off()
