set.seed(1)
# Load the required packages
library(tidyverse)
library(latex2exp)
library(GGally)
#specify options for saving the plots to files
single_width = 9
double_width = 14
single_height = 8.2
double_height = 12.8
res = 600

# Load results of the sensitivity analysis
sens_results = readr::read_csv(file = "sensitivity-analysis-results-main.csv")

# Path to save results
path = "Figures/main/"

# Histogram for R_h.
sens_results %>%
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
ggsave(filename = paste0(path, "colo_results_ICA.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)


# Histogram for Spearman's rho.
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
ggsave(filename = paste0(path, "colo_results_sprho.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Compute Estimated interval of ignorance.
range(sens_results$ICA)
range(sens_results$sp_rho)

# Proportions of dependent censoring of TTP by OS in both treatment groups.
mean(sens_results$prop_always + sens_results$prop_harmed)
mean(sens_results$prop_always + sens_results$prop_protected)

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


