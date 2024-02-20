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

# Load results of the sensitivity analyses. We do not use the main sensitivity
# analysis because the same results are available in
# sensitivity-analysis-results-relaxed.rds (together with sensitivity analyses
# under different assumptions).
sens_results_relaxed = read_rds("sensitivity-analysis-results-relaxed.rds")
sens_results_relaxed_tbl = sens_results_relaxed %>%
  rowwise() %>%
  reframe(sens_results) %>%
  select(-ranges, -mutinfo_estimator)
sens_results_main = sens_results_relaxed %>%
  filter(copula_family == "gaussian", range_class == "Main Assumptions", cond_ind)
sens_results_main_tbl = sens_results_relaxed_tbl %>%
  filter(copula_family == "gaussian", range_class == "Main Assumptions", cond_ind)

best_fitted_model = read_rds("best-fitted-model.rds")

# Path to save results
path_main = "Figures/main/"
path_relaxed = "Figures/relaxed/"
path_other_copulas = "Figures/other-copulas/"



# Histograms --------------------------------------------------------------

# Histogram for ICA = R_h excluding patients that die before progressing under
# both treatments; main analysis
sens_results_relaxed %>%
  filter(copula_family == "")
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
ggsave(filename = paste0(path, "Rh-subset.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Histogram for ICA = Spearman's rho excluding patients that die before
# progressing under both treatments.



# Histogram for ICA = Spearman's rho for the full population.
sens_results %>%
  ggplot(aes(x = sp_rho)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$\\rho_s$")) +
  geom_histogram(
    mapping = aes(y = after_stat(density)),
    fill = "gray",
    color = "black",
    binwidth = 0.025,
    boundary = 1
  ) +
  scale_y_continuous(name = "Density") +
  theme_bw()
ggsave(filename = paste0(path, "sp-rho-full.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)


# Uncertainty Intervals ---------------------------------------------------

# Compute Estimated intervals of ignorance and uncertainty.
sink(file = "sensitivity-intervals-main")
set.seed(1)
sensitivity_intervals_Rh_subset = sensitivity_intervals_Dvine(
  fitted_model = best_fitted_model,
  sens_results = sens_results,
  B = 200
)
set.seed(1)
sensitivity_intervals_sprho_full = sensitivity_intervals_Dvine(
  fitted_model = best_fitted_model,
  sens_results = sens_results,
  measure = "sp_rho",
  B = 200
)
# Dependent Censoring -----------------------------------------------------

# Proportions of dependent censoring of TTP by OS in both treatment groups.
mean(sens_results$prop_never + sens_results$prop_harmed)
mean(sens_results$prop_never + sens_results$prop_protected)

# Principal Strata of Diseases Status -------------------------------------











# Additional Results ------------------------------------------------------

# Additional exploration of how the assumptions translate to some easy to
# interpret quantities. We first look at survival classification probabilities
# in all simulated scenarios.
sens_results %>%
  pivot_longer(cols = starts_with("prop"),
               names_to = "type_prop",
               values_to = "Proportion") %>%
  mutate(type_prop = fct_recode(type_prop,
                                "Proportion Always" = "prop_always",
                                "Proportion Harmed" = "prop_harmed",
                                "Proportion Never" = "prop_never",
                                "Proportion Protected" = "prop_protected")) %>%
  ggplot(aes(x = Proportion)) +
  geom_histogram(fill = "gray", color = "black") +
  facet_wrap("type_prop", nrow = 2, ncol = 2) +
  scale_y_continuous(name = "Count") +
  theme_bw()
ggsave(filename = paste0(path, "survival_classification.png"),
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

sens_results %>%
  pivot_longer(cols = 4:7,
               names_to = "type_corr",
               values_to = "value") %>%
  mutate(type_corr = fct_recode(type_corr,
                                "S_0 and S_1" = "sp_rho_s0s1",
                                "S_0 and T_1" = "sp_rho_s0t1",
                                "S_1 and T_0" = "sp_rho_t0s1",
                                "T_0 and T_1" = "sp_rho_t0t1")) %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill = "gray", color = "black") +
  facet_wrap("type_corr", nrow = 2, ncol = 2) +
  scale_y_continuous(name = "Count") +
  scale_x_continuous(name = TeX("$\\rho_s$"), limits = c(0, 1)) +
  theme_bw()
ggsave(filename = paste0(path, "pairwise_sp_rho.png"),
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

sens_results %>%
  rename("S_0 and S_1" = "sp_rho_s0s1",
         "S_0 and T_1" = "sp_rho_s0t1",
         "S_1 and T_0" = "sp_rho_t0s1",
         "T_0 and T_1" = "sp_rho_t0t1") %>%
  ggpairs(columns = 4:7) +
  scale_x_continuous(name = TeX("$\\rho_s$")) +
  scale_y_continuous(name = TeX("$\\rho_s$")) +
  theme_bw()
ggsave(filename = paste0(path, "ggpairs_sp_rho.png"),
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

# Compute the range of the OR's for the principal strate of diseased status.
sens_results %>%
  mutate(OR = (prop_always * prop_never) / (prop_harmed * prop_protected)) %>%
  summarise(min_OR = min(OR), max_OR = max(OR))

