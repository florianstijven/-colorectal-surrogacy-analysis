# Setup -------------------------------------------------------------------
set.seed(1)
# Load the required packages
library(tidyverse)
library(latex2exp)
library(GGally)
library(Surrogate)
#specify options for saving the plots to files
single_width = 9
double_width = 14
single_height = 8.2
double_height = 12.8
res = 600

# The results of the main analysis are contained in two types of files: (i)
# results of the MC sensitivity analysis and (ii) uncertainty intervals. We ran
# two separate sensitivity analysis: one with the SICC as the ICA and one with
# Spearman's rho as the ICA.
sens_results_main_sicc = readRDS("results/sensitivity-analysis-results-main-sicc.rds")
sens_results_main_sprho = readRDS("results/sensitivity-analysis-results-main-sprho.rds")
# Sensitivity intervals where the ICA is defined in the population excluding
# patients where S_0 = T_0 AND S_1 = T_1.
sensitivity_intervals_sicc_subset = readRDS("results/sensitivity-intervals-Rh-subset.rds")
sensitivity_intervals_sprho_subset = readRDS("results/sensitivity-intervals-sprho-subset.rds")
# Sensitivity intervals where the ICA is defined in the full population. The
# SICC is not defined in this population, so we only use Spearman's rho here.
sensitivity_intervals_sprho_full = readRDS("results/sensitivity-intervals-sprho-full.rds")

# The results of the analyses under relaxed assumptions are contained in a
# single file. This file contains a tibble with one sensitivity analysis on each
# row.
sens_results_relaxed = readRDS("results/sensitivity-analysis-results-relaxed.rds")
# For plotting, the original tibble with one sensitivity analysis per row should
# be converted to a tibble with one MC replication per row.
sens_results_relaxed_tbl = sens_results_relaxed %>%
  rowwise() %>%
  reframe(sens_results) %>%
  select(-ranges, -mutinfo_estimator)

best_fitted_model = read_rds("results/best-fitted-model.rds")

# Path to save results that appear in the paper or its appendices.
save_to_main = "paper-figures-tables/main-text/"
save_to_appendix = "paper-figures-tables/appendix/"
# Other results are saved into a separate directory.
path_main = "additional-figures-tables/main-analysis/"
path_relaxed = "additional-figures-tables/relaxed-analysis/"

# Histograms --------------------------------------------------------------

# Histogram for ICA = R_h excluding patients that die before progressing under
# both treatments; main analysis.
sens_results_main_sicc %>%
ggplot(aes(x = ICA)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$R_h^2$")) +
  geom_histogram(
    mapping = aes(y = after_stat(density)),
    fill = "gray",
    color = "black",
    binwidth = 0.025,
    boundary = 1,
  ) +
  scale_y_continuous(name = "Density") +
  theme_bw()
ggsave(filename = paste0(save_to_main, "Rh-subset.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Histogram for ICA = R_h excluding patients that die before progressing under
# both treatments; relaxed analyses.
sens_results_relaxed_tbl %>%
  filter(ICA_type == "SICC") %>%
  ggplot(aes(x = ICA)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$R_h^2$")) +
  geom_histogram(
    mapping = aes(y = after_stat(density)),
    fill = "gray",
    color = "black",
    binwidth = 0.025,
    boundary = 1,
  ) +
  scale_y_continuous(name = "Density") +
  geom_histogram(data = tibble(ICA=c(-1, 2)), color = 'white') +
  facet_grid(copula_family~range_class) +
  theme_bw()
ggsave(filename = paste0(path_relaxed, "Rh-subset.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Histogram for ICA = Spearman's rho excluding patients that die before
# progressing under both treatments; main analysis.
sens_results_main_sprho %>%
  ggplot(aes(x = ICA)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$\\rho_{sp}$")) +
  geom_histogram(
    mapping = aes(y = after_stat(density)),
    fill = "gray",
    color = "black",
    binwidth = 0.025,
    boundary = 1,
  ) +
  scale_y_continuous(name = "Density") +
  theme_bw()
ggsave(filename = paste0(path_main, "sprho-subset.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Histogram for ICA = Spearman's rho excluding patients that die before
# progressing under both treatments; relaxed analyses.
sens_results_relaxed_tbl %>%
  filter(ICA_type == "Spearman's correlation") %>%
  ggplot(aes(x = ICA)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$\\rho_{sp}$")) +
  geom_histogram(
    mapping = aes(y = after_stat(density)),
    fill = "gray",
    color = "black",
    binwidth = 0.025,
    boundary = 1,
  ) +
  scale_y_continuous(name = "Density") +
  geom_histogram(data = tibble(ICA=c(-1, 2)), color = 'white') +
  facet_grid(copula_family~range_class) +
  theme_bw()
ggsave(filename = paste0(path_relaxed, "sprho-subset.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Histogram for ICA = Spearman's rho for the full population; main analysis.
sens_results_main_sicc %>%
  ggplot(aes(x = sp_rho)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$\\rho_{sp}$")) +
  geom_histogram(
    mapping = aes(y = after_stat(density)),
    fill = "gray",
    color = "black",
    binwidth = 0.025,
    boundary = 1,
  ) +
  scale_y_continuous(name = "Density") +
  theme_bw()
ggsave(filename = paste0(path_main, "sprho-full.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Histogram for ICA = Spearman's rho for the full population; relaxed analyses.
# The results for sp_rho in the full population are repeated twice, once for
# each value of ICA_type. We therefore select one ICA_type.
sens_results_relaxed_tbl %>%
  filter(ICA_type == "Spearman's correlation") %>%
  ggplot(aes(x = sp_rho)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$\\rho_{sp}$")) +
  geom_histogram(
    mapping = aes(y = after_stat(density)),
    fill = "gray",
    color = "black",
    binwidth = 0.025,
    boundary = 1,
  ) +
  scale_y_continuous(name = "Density") +
  geom_histogram(data = tibble(sp_rho=c(-1, 2)), color = 'white') +
  facet_grid(copula_family~range_class) +
  theme_bw()
ggsave(filename = paste0(path_relaxed, "sprho-full.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)


# Uncertainty Intervals ---------------------------------------------------

# To prevent having to run this file mutliple times, the sensitivity intervals
# are printed to a .txt file.
sink(file = paste0(save_to_main, "sensitivity-intervals.txt"))
cat(paste(rep("=", 80), collapse = '')); cat("\n\n")
cat("ICA  = SICC in the subpopulation of patients where patients with S_0 = T_0 AND S_1 = T_1 are excluded.\n\n")
print(sensitivity_intervals_sicc_subset)
cat("\n"); cat(paste(rep("=", 80), collapse = '')); cat("\n\n")

cat("ICA  = Spearman's rho in the subpopulation of patients where patients with S_0 = T_0 AND S_1 = T_1 are excluded.\n\n")
print(sensitivity_intervals_sprho_subset)
cat("\n"); cat(paste(rep("=", 80), collapse = '')); cat("\n\n")

cat("ICA  = Spearman's rho in the full population.\n\n")
print(sensitivity_intervals_sprho_full)
cat("\n"); cat(paste(rep("=", 80), collapse = '')); cat("\n\n")
sink()

# For the sensitive intervals under relaxed assumptions, we add the sensitivity
# intervals directly to the original tibble.
interval_paste = function(interval, digits = 3) {
  paste0("(", round(interval[1], digits), ", ", round(interval[2], digits), ")")
}
sink(file = paste0(save_to_appendix, "sensitivity-intervals-relaxed.txt"))
sens_results_relaxed %>%
  ungroup() %>%
  rowwise(c("range_class", "copula_family", "ICA_type")) %>%
  summarise(
    "Est. II, ICA in subset" = interval_paste(sens_interval_ICA_subset$est_interval_of_ignorance, 2),
    "UI (pointwise), ICA in subset" = interval_paste(
      sens_interval_ICA_subset$interval_of_uncertainty_pointwise_coverage,
      2
    ),
    "Est. II, ICA = sp_rho in full pop." = interval_paste(sens_interval_sprho_full$est_interval_of_ignorance, 2),
    "UI (pointwise), ICA = sp_rho in full pop." = interval_paste(
      sens_interval_sprho_full$interval_of_uncertainty_pointwise_coverage,
      2
    ),
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = "ICA_type", values_from = 4:7) %>%
  select(c(1, 2, 3, 5, 4, 6, 7, 9)) %>%
  knitr::kable(format = "latex") %>%
  print()
sink()

# Dependent Censoring -----------------------------------------------------

# Proportions of dependent censoring of TTP by OS in both treatment groups:
# P(S_0 = T_0) = 0.333 and P(S_1 = T_1) = 0.304.
sink(file = paste0(save_to_appendix, "proportion-dependent-censoring.txt"))
cat("Estimated proportion dependent censoring for Z = 0: P(S_0 = T_0) = ")
cat(mean(sens_results_main_sicc$prop_never + sens_results_main_sicc$prop_harmed))
cat("\n")
cat("Estimated proportion dependent censoring for Z = 1: P(S_1 = T_1) = ")
cat(mean(sens_results_main_sicc$prop_never + sens_results_main_sicc$prop_protected))
sink()

# Additional Results ------------------------------------------------------

# Additional exploration of how the assumptions translate to some easy to
# interpret quantities. We first look at survival classification probabilities
# in all simulated scenarios.
sink(file = paste0(save_to_appendix, "ranges-population-strata.txt"))
sens_results_main_sicc %>%
  pivot_longer(cols = starts_with("prop"),
               names_to = "type_prop",
               values_to = "Proportion") %>%
  group_by(type_prop) %>%
  summarise(min_prop = min(Proportion), max_prop = max(Proportion)) %>%
  print()
sink()

# Model-based Spearman's rho between S_0 and T_0, and between S_1 and T_1: 0.673
# and 0.676, respectively.
mean(sens_results_main_sicc$sp_rho_t0s0); mean(sens_results_main_sicc$sp_rho_s1t1)
sp_rho_observable = 0.5 * (
  mean(sens_results_main_sicc$sp_rho_t0s0) + mean(sens_results_main_sicc$sp_rho_s1t1)
)

sens_results_main_sicc %>%
  rename(
    "S_0 and S_1" = "sp_rho_s0s1",
    "S_0 and T_1" = "sp_rho_s0t1",
    "S_1 and T_0" = "sp_rho_t0s1",
    "T_0 and T_1" = "sp_rho_t0t1"
  ) %>%
  ggpairs(
    columns = 4:7,
    diag = list(continuous = wrap("densityDiag", rescale = TRUE)),
    upper = list(continuous = "blank"),
    lower = list(
      continuous = function(data, mapping, ...) {
        fn = wrap("points")
        p = fn(data, mapping, ...)
        p + geom_vline(
          xintercept = sp_rho_observable,
          color = "red",
          alpha = 0.50
        ) +
          geom_hline(yintercept = sp_rho_observable,
                     color = "red",
                     alpha = 0.50)
      }
    )
  ) +
  scale_x_continuous(
    name = TeX("$\\rho_s$"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  xlab(TeX("$\\rho_{sp}(X, Y)$")) +
  ylab(TeX("$\\rho_{sp}(X, Y)$")) +
  scale_y_continuous(
    name = TeX("$\\rho_s$"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  theme_bw()
ggsave(filename = paste0(save_to_appendix, "ggpairs_sp_rho.png"),
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

# Compute the range of the OR's for the principal strata of diseased status.
sink(file = paste0(path_relaxed, "OR-diseased-status.txt"))
sens_results_relaxed_tbl %>%
  mutate(OR = (prop_always * prop_never) / (prop_harmed * prop_protected)) %>%
  group_by(range_class, copula_family, ICA_type) %>%
  summarize("Mininimum OR" = min(OR), "Maximum OR" = max(OR)) %>%
  ungroup()
sink()

