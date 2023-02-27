set.seed(1)
# Load the required packages
library(tidyverse)
library(latex2exp)
library(GGally)
#set up parameters for saving plots
path = "Code - Combining Conditions/"
width = 13.2
height = 9.9
pointsize = 10
res = 500

# Load results of the sensitivity analysis
load(file = "Code - Combining Conditions/sensitivity_analysis_results_100k.RData")

# Imposing one restriction at a time is often an insufficient representation of
# the available subject-matter knowledge. Therefore, it can be of interest to
# impose different types of assumptions simultaneously.

# As an illustration, we look at 4 possible combinations of the following
# assumptions:
# * weaker cross-associations and stochastic monotonicity 
# * odds ratio of diseased in (0, 4) or larger than 3
# * within associations in (0.5, 0.975) or (0.5, 0.85)(in terms of Spearman's rho)

# These restrictions are imposed by filtering the data frame from the
# sensitivity analysis.

# Construct an artificial data set to ease further plots.
sens_results_helper = sens_results_no_cond_ind %>%
  filter(pmin(abs(sp_s0t0), abs(sp_s1t1)) > pmax(abs(sp_s0t1), abs(sp_s1t0))) %>%
  filter(sp_s0s1 > 0,
         sp_s0t0 > 0,
         sp_s0t1 > 0,
         sp_s1t0 > 0,
         sp_s1t1 > 0,
         sp_t0t1 > 0) %>%
  filter((prop_always * prop_never) / (prop_harmed * prop_protected) < 4) %>%
  filter((prop_always * prop_never) / (prop_harmed * prop_protected) > 0) %>%
  filter(sp_s0s1 > 0.5, sp_t0t1 > 0.5) %>%
  filter(sp_s0s1 < 0.975, sp_t0t1 < 0.975) %>%
  mutate(
    odds_ratio = factor(
      x = 1,
      levels = 1:2,
      labels = c("(0, 4)", "> 4")
    ),
    within_associations = factor(
      x = 1,
      levels = 1:2,
      labels = c("(0.5, 0.975)", "(0.5, 0.85)")
    )
  )


sens_results_helper = sens_results_helper %>%
  bind_rows(
    sens_results_no_cond_ind %>% 
      filter(pmin(abs(sp_s0t0), abs(sp_s1t1)) > pmax(abs(sp_s0t1), abs(sp_s1t0))) %>%
      filter(
        sp_s0s1 > 0,
        sp_s0t0 > 0,
        sp_s0t1 > 0,
        sp_s1t0 > 0,
        sp_s1t1 > 0,
        sp_t0t1 > 0
      ) %>%
      filter((prop_always * prop_never) / (prop_harmed * prop_protected) > 4) %>%
      filter(sp_s0s1 > 0.5, sp_t0t1 > 0.5) %>%
      filter(sp_s0s1 < 0.975, sp_t0t1 < 0.975) %>%
      mutate(
        odds_ratio = factor(
          x = 2,
          levels = 1:2,
          labels = c("(0, 4)", "> 4")
        ),
        within_associations = factor(
          x = 1,
          levels = 1:2,
          labels = c("(0.5, 0.975)", "(0.5, 0.85)")
        )
      )
  )

sens_results_helper = sens_results_helper %>%
  bind_rows(
    sens_results_no_cond_ind %>% 
      filter(pmin(abs(sp_s0t0), abs(sp_s1t1)) > pmax(abs(sp_s0t1), abs(sp_s1t0))) %>%
      filter(
        sp_s0s1 > 0,
        sp_s0t0 > 0,
        sp_s0t1 > 0,
        sp_s1t0 > 0,
        sp_s1t1 > 0,
        sp_t0t1 > 0
      ) %>%
      filter((prop_always * prop_never) / (prop_harmed * prop_protected) > 4) %>%
      filter(sp_s0s1 > 0.5, sp_t0t1 > 0.5) %>%
      filter(sp_s0s1 < 0.85, sp_t0t1 < 0.85) %>%
      mutate(
        odds_ratio = factor(
          x = 2,
          levels = 1:2,
          labels = c("(0, 4)", "> 4")
        ),
        within_associations = factor(
          x = 2,
          levels = 1:2,
          labels = c("(0.5, 0.975)", "(0.5, 0.85)")
        )
      )
  )

sens_results_helper = sens_results_helper %>%
  bind_rows(
    sens_results_no_cond_ind %>% 
      filter(pmin(abs(sp_s0t0), abs(sp_s1t1)) > pmax(abs(sp_s0t1), abs(sp_s1t0))) %>%
      filter(
        sp_s0s1 > 0,
        sp_s0t0 > 0,
        sp_s0t1 > 0,
        sp_s1t0 > 0,
        sp_s1t1 > 0,
        sp_t0t1 > 0
      ) %>%
      filter((prop_always * prop_never) / (prop_harmed * prop_protected) < 4) %>%
      filter((prop_always * prop_never) / (prop_harmed * prop_protected) > 0) %>%
      filter(sp_s0s1 > 0.5, sp_t0t1 > 0.5) %>%
      filter(sp_s0s1 < 0.975, sp_t0t1 < 0.975) %>%
      mutate(
        odds_ratio = factor(
          x = 1,
          levels = 1:2,
          labels = c("(0, 4)", "> 4")
        ),
        within_associations = factor(
          x = 2,
          levels = 1:2,
          labels = c("(0.5, 0.975)", "(0.5, 0.85)")
        )
      )
  )


# plot the resulting histograms
png(
  filename = paste0(path, "colo_combining_conditions_sprho.png"),
  width = width,
  height = height,
  units = "cm",
  pointsize = pointsize,
  res = res
)
sens_results_helper %>%
  ggplot(aes(x = sp_rho)) + 
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$\\rho_s$")) +
  geom_histogram(
    fill = "gray",
    color = "black",
    binwidth = 0.025,
    boundary = 1
  ) +
  theme_bw() +
  facet_grid(rows = vars(within_associations), 
             cols = vars(odds_ratio)
             )
dev.off()
