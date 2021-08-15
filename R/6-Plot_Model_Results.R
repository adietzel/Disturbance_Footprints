### Plot model results ################################################
### Author: Andreas Dietzel, andreas.dietzel@my.jcu.edu.au ############

rm(list = ls())

# Load libraries
library(tidyverse)

# Custom theme for plots
theme_simple <- theme_classic() +
  theme(text       = element_text(color = "black"),
        strip.text = element_text(color = "black"),
        axis.text  = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        line       = element_line(color = "black"),
        plot.background   = element_rect(fill = "white", color = "transparent"),
        panel.background  = element_rect(fill = "white", color = "black"),
        strip.background  = element_rect(fill = "white", color = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.background = element_rect(fill = "white", color = "transparent"),
        legend.key        = element_rect(fill = "white", color = "transparent"))

# Load data
Reefs <- read.csv("data/Centroids.csv", strip.white = T) # reef locations

# Model output files
Impact.small <- read.csv("RData/Impact_results_Small.csv") %>%
  set_names(c("Autocorrelation", "Species", "Iter", "Impact")) %>%
  mutate(Magnitude = "Small (5% of reefs disturbed)")
Impact.med <- read.csv("RData/Impact_results_Medium.csv") %>%
  set_names(c("Autocorrelation", "Species", "Iter", "Impact")) %>%
  mutate(Magnitude = "Medium (25% of reefs disturbed)")
Impact.large <- read.csv("RData/Impact_results_Large.csv") %>%
  set_names(c("Autocorrelation", "Species", "Iter", "Impact")) %>%
  mutate(Magnitude = "Large (50% of reefs disturbed)")

# Row-bind impact results and rename species and degree of autocorrelation
Impact <- bind_rows(Impact.small, Impact.med, Impact.large) %>%
  mutate(Autocorrelation = ifelse(Autocorrelation == 1, "no autocorrelation (random)",
                                  ifelse(Autocorrelation == 2, "low autocorrelation",
                                         "high autocorrelation")),
         Species = ifelse(Species == 1, "Global disperser",
                          ifelse(Species == 2, "Long-distance spawner",
                                 ifelse(Species == 3, "Short-distance spawner", "Brooder"))))

# Calculate relative impact by simulation run (Iter)
Impact.relative <- Impact %>%
  group_by(Autocorrelation, Magnitude, Iter) %>%
  mutate(Impact.relative = Impact/min(Impact))

# Summarise mean and standard deviation of relative impact for each species, degree of autocorrelation and magnitude
Impact.summary <- as.data.frame(Impact.relative) %>%
  group_by(Autocorrelation, Species, Magnitude) %>%
  summarise(Mean = mean(Impact.relative),
            Lower = quantile(Impact.relative, probs = .05),
            Upper = quantile(Impact.relative, probs = .95),
            .groups = "keep")

# Set factor levels for species
Impact.summary$Species <- factor(Impact.summary$Species,
                                 levels = c("Global disperser", "Long-distance spawner",
                                            "Short-distance spawner", "Brooder"))

# Set factor levels for magnitude
Impact.summary$Magnitude <- factor(
  Impact.summary$Magnitude,
  levels = c("Small (5% of reefs disturbed)", "Medium (25% of reefs disturbed)",
             "Large (50% of reefs disturbed)"))

# Set factor levels for autocorrelation
Impact.summary$Autocorrelation <- factor(Impact.summary$Autocorrelation,
                                         levels = c("no autocorrelation (random)",
                                                    "low autocorrelation",
                                                    "high autocorrelation"))

# Plot relative impact
(Impact.boxplot <- ggplot(Impact.summary, aes(x = Species, y = Mean, colour = Autocorrelation)) +
    geom_pointrange(aes(ymin = Lower, ymax = Upper),
                    position = position_dodge(width = .7)) +
    facet_wrap(~ Magnitude, ncol = 1) +
    theme_simple +
    geom_vline(xintercept = c(1.5,2.5,3.5), colour = "lightgrey", linetype = 3) +
    labs(x = "", colour = "",
         y = expression(paste("Relative impact ( ",
                              integral(Delta),"N"["t  (L,S,B)"], " dt ", "   /  ",
                              integral(Delta),"N"["t  (G)"], " dt ", " )"))) +
    scale_colour_manual(values = c("#EB563A", "#44A9CC", "#FFB600")) +
    scale_x_discrete(labels = function(x){sub("\\s", "\n", x)}) +
    theme(legend.position = c(.28,.95),
          legend.background = element_blank()))

# Save plot
ggsave("figures/Figure5_ImpactPlot.jpg", width = 4.5, height = 7, dpi = 600)
