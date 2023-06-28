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
sens_results = readr::read_csv(file = "VSC/results/sensitivity-analysis-results-relaxed.csv")

# Path to save results
path = "Figures/relaxed/"

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
