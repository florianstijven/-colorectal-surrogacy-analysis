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

# Load results of the sensitivity analysis. We combine the main results (w/
# Gaussian copulas) and the results where the unidentifiable copulas are of a
# different type.
sens_results = bind_rows(
  readr::read_csv(file = "VSC/results/sensitivity-analysis-results-other-copulas.csv"),
  readr::read_csv(file = "VSC/results/sensitivity-analysis-results-main.csv") %>%
    mutate(copula = "Gaussian")
)

# Path to save results
path = "Figures/other-copulas/"

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
  facet_wrap("copula") +
  theme_bw()
ggsave(filename = paste0(path, "colo_results_ICA.png"),
       device = "png",
       width = double_width,
       height = double_height,
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
  facet_wrap("copula") +
  theme_bw()
ggsave(filename = paste0(path, "colo_results_sprho.png"),
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Compute Estimated interval of ignorance.
sens_results %>%
  group_by(copula) %>%
  summarise(min_ICA = min(ICA), max_ICA = max(ICA),
            min_sp_rho = min(sp_rho), max_sp_rho = max(sp_rho))

# Plot to assess plausibility of the results.
sens_results %>%
  ggplot(aes(x = sp_rho_s0t1, y = ICA)) +
  geom_point() +
  facet_grid("copula")

