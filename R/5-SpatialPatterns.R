### Generate spatial disturbance patterns #############################
### Author: Andreas Dietzel, andreas.dietzel@my.jcu.edu.au ############

rm(list = ls())

# Load libraries
library(sp)
library(gstat)
library(sf)
library(ncf)
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
Reefs <- read.csv("data/Centroids.csv", strip.white = T) %>%
  filter(Lon > 0) # reef locations

# Create scenario overview
scenarios <- data.frame(
  Extent = c(rep("small", 3), rep("medium", 3), rep("large", 3)),
  Autocor = c(rep(c("no","low","high"), 3)),
  Range = c(.05, .1, 50, .05, .25, 50, .05, .5, 500),
  Magn = c(rep(round(.05*nrow(Reefs)), 3),
           rep(round(.25*nrow(Reefs)), 3),
           rep(round(.5*nrow(Reefs)), 3))
)

# Create empty data frames for storing disturbance patterns
df <- bind_cols(Reefs, data.frame(matrix(ncol = 100*nrow(scenarios), nrow = nrow(Reefs))))
df2 <- bind_cols(Reefs, data.frame(matrix(ncol = 100*nrow(scenarios), nrow = nrow(Reefs))))

# Define for each scenario an exponential variogram model to predict for each reef whether it gets disturbed or not
for (i in 1:nrow(scenarios)) {

  # Define variogram model for each scenario
  vario <- vgm(psill = 1, range = scenarios[i, "Range"], model = 'Exp')

  zDummy <- gstat(formula = z~1, locations = ~Lon+Lat, dummy = TRUE,
                  beta = 1, model = vario, nmax = 20)

  # Use variogram model to predict 100 disturbance patterns for each scenario
  df2[, c(((i-1) * 100 + 3):(i * 100 + 2))] <-
    predict(zDummy, newdata = Reefs, nsim = 100) %>%
    dplyr::select(-Lon,-Lat)

  # Name columns programmatically for each scenario and simulation run
  colnames(df2)[c(((i-1) * 100 + 3):(i * 100 + 2))] <-
    paste(scenarios[i, "Extent"], scenarios[i, "Autocor"], 1:100, sep = "_")

  # Convert continuous variogram prediction variable into binary (disturbed vs not disturbed)
  df[, c(((i-1) * 100 + 3):(i * 100 + 2))] <-
    df2[, c(((i-1) * 100 + 3):(i * 100 + 2))] %>%
    mutate(across(where(is.numeric),
                  ~ ifelse(.x < stats::quantile(.x, probs = scenarios[i, "Magn"]/nrow(Reefs)), 1, 0)))

  # Rename columns
  colnames(df)[c(((i-1) * 100 + 3):(i * 100 + 2))] <-
    paste(scenarios[i, "Extent"], scenarios[i, "Autocor"], 1:100, sep = "_")
}

# Save file with simulated spatial patterns to use in Matlab
# write.csv(df, "figures/SpatialPatterns.csv")

# Create example maps for each combo of magnitude x autocorrelation
Maps.df <- df %>%
  select(Lat, Lon, ends_with("_10")) %>%
  pivot_longer(-c(1,2), names_to = "Event", values_to = "Value") %>%
  separate(Event, c("Magnitude", "Autocorrelation", "Iter"), "_") %>%
  mutate(Value2 = as.factor(ifelse(Value == 0, "a", Autocorrelation))) #

# Reorder factor levels
Maps.df$Magnitude <- factor(Maps.df$Magnitude, levels = c("small","medium","large"))
Maps.df$Autocorrelation <- factor(Maps.df$Autocorrelation, levels = c("no","low","high"))

# Plot the maps
(Maps.plot <-
    ggplot() +
    geom_point(data = Maps.df,
               aes(x = Lon, y = Lat, colour = Value2), size = .0001) +
    coord_map() +
    theme_void() + theme(legend.position = "") +
    facet_grid(Magnitude ~ Autocorrelation) +

    scale_colour_manual(values = c("gray93", "#FFB600", "#44A9CC", "#EB563A")))

# Generate examplary correlograms for figure 4

a2 <- df %>%
  select(Lat, Lon, ends_with("_10")) %>%
  pivot_longer(-c(1,2), names_to = "Event", values_to = "Value") %>%
  separate(Event, c("Magnitude", "Autocorrelation", "Iter"), "_")

# Reorder factor levels
a2$Magnitude <- factor(a2$Magnitude, levels = c("small","medium","large"))
a2$Autocorrelation <- factor(a2$Autocorrelation, levels = c("high","low","no"))

Corr.df <- data.frame()

for (i in 1:nrow(scenarios)) {


  xyz <- a2 %>% filter(Magnitude == scenarios[i, "Extent"]) %>%
    filter(Autocorrelation == scenarios[i, "Autocor"])

  # Generate spline cross-correlogram
  Corr <- spline.correlog(xyz$Lon, xyz$Lat, xyz$Value,
                          npoints = 100, resamp = 100, latlon = T)

  # Store output of cross-correlogram
  Corr.df.temp <- data.frame(distance = t(Corr$real$predicted$x),
                             correlation = t(Corr$real$predicted$y),
                             lower = t(Corr$boot$boot.summary$predicted$y)[,2],
                             upper = t(Corr$boot$boot.summary$predicted$y)[,10],
                             Autocorrelation = scenarios[i, "Autocor"],
                             Magnitude = scenarios[i, "Extent"])

  # Bind rows
  Corr.df <- bind_rows(Corr.df, Corr.df.temp)
}

# Reorder factor levels
Corr.df$Magnitude <- factor(Corr.df$Magnitude, levels = c("small","medium","large"))
Corr.df$Autocorrelation <- factor(Corr.df$Autocorrelation, levels = c("no","low","high"))

# Plot the correlograms
(Corr.plot.large <- Corr.df %>%
    ggplot(aes(x = distance, y = correlation)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Autocorrelation),
                alpha = .3, colour = NA) +
    geom_line(aes(colour = Autocorrelation)) +
    coord_cartesian(xlim = c(0, 500), ylim = c(-.3, 1)) +
    labs(linetype = "", x = "Distance (km)", y = "Correlation",
         colour = "") +
    theme_simple +
    facet_wrap(~ Magnitude, ncol = 1) +
    scale_colour_manual(values = c("#EB563A", "#44A9CC", "#FFB600")) +
    scale_fill_manual(values = c("#EB563A", "#44A9CC", "#FFB600")) +
    theme(legend.position = ""))

# Bring all plots together
cowplot::plot_grid(Impact.boxplot, Maps.plot, Corr.plot.large,
                   rel_widths = c(.7, 1, .5), labels = "auto", nrow = 1)

ggsave("figures/ModelOutputFigure4.pdf", width = 14, height = 8)

# Save output files
save(list = ls(), file = "RData/AllFiles.RData")

