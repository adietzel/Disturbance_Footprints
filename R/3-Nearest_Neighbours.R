### Stats on spatial autocorrelation nearest neighbours ###############
### Author: Andreas Dietzel, andreas.dietzel@my.jcu.edu.au ############

### Calculate two clustering metrics:
### a) distance to nearest undisturbed neighbour
### b) proportion of undisturbed sites within 100km of disturbed site

rm(list = ls())

# LOAD LIBRARIES AND DATA

library(tidyverse)
library(cowplot)
library(fields)
library(stringr)
library(tidybayes)
library(modelr)
library(emmeans)
library(brms)

load("RData/0-Events.RData")
load("RData/0-Auxiliary.RData")
load("RData/2-SplineCor.RData")
load("RData/3-NN.RData")
theme_set(theme_simple)

# CALCULATE DISTANCE METRICS FOR EACH EVENT (SEVERITY THRESHOLD > 2)

a <- lapply(seq_along(unique(Events.Agg.df$Event)), function(i) {

  data = Events.Agg.df[Events.Agg.df$Event == unique(Events.Agg.df$Event)[i],] %>%
    mutate(ID1 = as.character(row_number()))

  coords = as.matrix(data.frame(Lon = data$Lon, Lat = data$Lat))

  distance.matrix <- as.data.frame(rdist.earth(coords, coords,
                                               miles = F, R = NULL))

  dist.mat.long <- distance.matrix %>%
    mutate(ID1 = as.character(1 : length(data$Lon))) %>%
    tidyr::gather(., ID2, distance, c(1 : length(data$Severe345))) %>%
    mutate(ID2 = substr(ID2, 2, 5),
           Self = as.numeric(ID1) - as.numeric(ID2)) %>%
    filter(!Self == 0) %>%
    left_join(., data[, c("ID1", "Severe345")], by = "ID1") %>%
    rename(Severe345A = Severe345) %>%
    left_join(., data[, c("ID1", "Severe345")], by = c("ID2" = "ID1")) %>%
    rename(Severe345B = Severe345)

  NN <- dist.mat.long %>% filter(Severe345B < 1) %>%
    filter(Severe345A > 0) %>% group_by(ID1) %>%
    summarise(MinDist = min(distance), .groups = "keep")

  PropN <- dist.mat.long %>% filter(distance < 100) %>%
    filter(Severe345A > 0) %>%
    group_by(ID1, Severe345B) %>%
    summarise(Freq = n(), .groups = "keep") %>% ungroup() %>%
    complete(ID1, Severe345B, fill = list(Freq = 0)) %>%
    group_by(ID1) %>%
    mutate(Prop = Freq/sum(Freq),
           FreqAll = sum(Freq)) %>%
    filter(Severe345B < 1) %>%
    dplyr::select(ID1, Prop, FreqAll, Freq) %>%
    left_join(., data[,c("Lon","Lat","ID1")], by="ID1")

  print(i)

  output = data.frame(Event = unique(Events.Agg.df$Event)[i],
                      MinDist = NN$MinDist,
                      PropN = PropN$Prop,
                      FreqAll = PropN$FreqAll,
                      Freq = PropN$Freq,
                      Lon = PropN$Lon,
                      Lat = PropN$Lat)
})

NN.emp <- bind_rows(a) %>%
  mutate(type = "Empirical",
         severe = ">2")

### CREATE NULL MODEL DISTRIBUTIONS
### Shuffle column Binary n times
### Record for each event the mean distance and proportion

sims = 100

b <- lapply(seq_along(unique(Events.Agg.df$Event)), function(i) {

  data = Events.Agg.df[Events.Agg.df$Event == unique(Events.Agg.df$Event)[i], ] %>%
    mutate(ID1 = as.character(row_number()))

  coords = as.matrix(data.frame(Lon = data$Lon, Lat = data$Lat))

  distance.matrix <- as.data.frame(rdist.earth(coords, coords,
                                               miles = F, R = NULL))

  dist.mat.long <- distance.matrix %>%
    mutate(ID1 = as.character(1 : length(data$Lon))) %>%
    tidyr::gather(., ID2, distance, c(1 : length(data$Severe345))) %>%
    mutate(ID2 = substr(ID2, 2, 5),
           Self = as.numeric(ID1) - as.numeric(ID2)) %>%
    filter(!Self == 0
           # & distance > 0
    ) %>%
    #### Remove distance > 0 because Yasi scores have same coordinates
    left_join(., data[, c("ID1", "Severe345")], by = "ID1") %>%
    rename(Severe345A = Severe345) %>%
    left_join(., data[, c("ID1", "Severe345")], by = c("ID2" = "ID1")) %>%
    rename(Severe345B = Severe345)

  MinDist = rep(NA, sims)
  PropN = rep(NA, sims)

  for (j in 1 : sims) {
    dist.mat.long$Severe345B <- sample(dist.mat.long$Severe345B, replace = FALSE)

    NN <- dist.mat.long %>% filter(Severe345B < 1) %>%
      filter(Severe345A > 0) %>% group_by(ID1) %>%
      summarise(MinDist = min(distance), .groups = "keep")

    MinDist <- c(MinDist, NN$MinDist)

    PropNN <- dist.mat.long %>% filter(distance < 100) %>%
      filter(Severe345A > 0) %>%
      group_by(ID1, Severe345B) %>%
      summarise(Freq = n(), .groups = "keep") %>%
      ungroup() %>%
      complete(ID1, Severe345B, fill = list(Freq = 0)) %>%
      group_by(ID1) %>%
      mutate(Prop = Freq / sum(Freq)) %>%
      filter(Severe345B < 1) %>%
      dplyr::select(ID1, Prop) %>%
      left_join(., data[,c("Lon", "Lat", "ID1")], by = "ID1")

    PropN <- c(PropN, PropNN$Prop)
  }

  print(i)

  Null = data.frame(Event = unique(Events.Agg.df$Event)[i],
                    MinDist = MinDist,
                    PropN = PropN)
})

NN.Null <- bind_rows(b) %>%
  filter(!is.na(MinDist)) %>%
  mutate(type = "Null", Lon = NA, Lat = NA,
         severe = ">2")

### MODEL FITS FOR NULL EXPECTATIONS

data.MinDist <- NN.Null %>%
  mutate(MinDist = ifelse(MinDist == 0, 0.01, MinDist))

model.MinDist.Null <- brm(log(MinDist) ~ Event, data = data.MinDist)

fits.MinDist.Null = NN.Null %>%
  data_grid(Event) %>%
  add_fitted_draws(model.MinDist.Null) %>%
  mutate(MinDist = exp(.value)) %>%
  group_by(Event) %>%
  summarise(MinDist.Null = median(MinDist), .groups = "keep") %>%
  mutate(severe = ">2")

data.PropN <- NN.Null %>%
  mutate(PropN = ifelse(PropN == 0, PropN + 0.001,
                        ifelse(PropN == 1, PropN - 0.001, PropN)))

model.PropN.Null <- brm(logit_scaled(PropN) ~ Event, data = data.PropN)

fits.PropN.Null = NN.Null %>%
  data_grid(Event) %>%
  add_fitted_draws(model.PropN.Null) %>%
  mutate(PropN = inv_logit_scaled(.value)) %>%
  group_by(Event) %>%
  summarise(PropN.Null = median(PropN), .groups = "keep") %>%
  mutate(severe = ">2")

### MODELS FITS FOR EMPIRICAL DATA

NN.emp2 <- NN.emp %>%
  mutate(PropN = ifelse(PropN == 0, PropN + 0.001,
                        ifelse(PropN == 1, PropN - 0.001, PropN)))

model.PropN.emp <- brm(bf(logit_scaled(PropN) ~ Event,
                          sigma ~ Event),
                       data = NN.emp2, family = "gaussian")

fits.PropN <- NN.emp2 %>%
  data_grid(Event) %>%
  add_fitted_draws(model.PropN.emp) %>%
  rename(PropN = .value) %>%
  mutate(PropN = inv_logit_scaled(PropN),
         severe = ">2")

model.MinDist.emp <- brm(log(MinDist) ~ Event, data = NN.emp2,
                         family = "gaussian",
                         control = list(adapt_delta = 0.995))

fits.NN <- NN.emp2 %>%
  data_grid(Event) %>%
  add_fitted_draws(model.MinDist.emp) %>%
  rename(NN = .value) %>%
  mutate(severe = ">2")

### PLOTTING PROPORTIONS AND MINIMUM DISTANCES TO PRODUCE FIGURE 2

(ModelFit.PropN.plot <- NN.emp %>%
    ggplot(aes(y = PropN, x = Event)) +
    geom_boxplot(outlier.shape = NA, alpha = .7, width = .5,
                 position = position_nudge(x = -0.2), fill = "lightgrey") +
    stat_pointinterval(aes(y = PropN), data = fits.PropN, .width = c(.66, .95),
                       position = position_nudge(x = 0.2), size = .7,
                       colour = "black") +
    geom_point(data = fits.PropN.Null,
               aes(x = Event, y = PropN.Null), colour = "red",
               position = position_nudge(x = 0.2)) +
    labs(x = "", y = "Proportion", colour = "Interval") +
    theme(axis.text.x = element_text(angle = 90, vjust = .5),
          legend.position = ""))

(ModelFit.MinDist.plot <- NN.emp %>%
    ggplot(aes(y = MinDist, x = Event)) +
    geom_boxplot(outlier.shape = NA, alpha = .7, width = .5,
                 position = position_nudge(x = -0.2), fill = "lightgrey") +
    stat_pointinterval(aes(y = exp(NN)), data = fits.NN, .width = c(.66, .95),
                       position = position_nudge(x = 0.2), size = .7,
                       colour = "black") +
    geom_point(data = fits.MinDist.Null,
               aes(x = Event, y = MinDist.Null), colour = "red",
               position = position_nudge(x = 0.2)) +
    scale_y_continuous(limits = c(0,350)) +
    guides(fill = F) +
    labs(x = "", y = "Distance (km)", colour = "Prediction\ninterval") +
    theme(axis.text.x = element_text(angle = 90, vjust = .5)))

### PLOT CORRELOGRAMS TO PRODUCE FIGURE 2

# # Dispersal kernels
# p_best <- read.csv("data/DispersalMike.csv", strip.white = T) %>%
#   mutate(Species = as.factor(1:1536)) %>%
#   filter(b > -1)
#
# library(plyr)
# coeflines <-
#   alply(as.matrix(p_best[,1:2]), 1, function(coef) {
#     stat_function(fun = function(x){exp(coef[2] * x) * coef[1] * (-5)},
#                   colour = "darkblue", alpha = .05, size = .1)
#   })
#
SplineCor.plot <- SC.AggReef %>%
  group_by(event) %>%
  filter(severe == ">2") %>%
  filter(distance < 0.5 * max(distance)) %>%
  ggplot(aes(x = distance, y = correlation, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(alpha = .3) +
  geom_line() +
  # coeflines +
  facet_wrap(~ event, ncol = 2) +
  coord_cartesian(xlim = c(0, 350), ylim = c(-.3, 1)) +
  labs(x = "Distance (km)", y = "Correlation")
#
# coeflines2 <-
#   alply(as.matrix(p_best[,1:2]), 1, function(coef) {
#     stat_function(fun = function(x){exp(coef[2] * x) * coef[1]},
#                   colour = "darkblue", alpha = .05, size = .1)
#   })
#
# DispKernels.plot <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
#   coeflines2 +
#   xlim(0,350) + theme_classic() +
#   coord_flip() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank()) +
#   labs(x = "", y = "kernel density")
#
# p2 <- insert_yaxis_grob(ModelFit.MinDist.plot,
#                         DispKernels.plot,
#                         grid::unit(.6, "null"), position = "right")

cowplot::plot_grid(
  SplineCor.plot, ncol = 1, labels = c("(a)",""),
  label_fontface = "italic", rel_heights = c(1.8,1),
  cowplot::plot_grid(ModelFit.PropN.plot, ModelFit.MinDist.plot,
                     nrow = 1,
                     labels = c("(b)","(c)"), label_fontface = "italic"))

ggsave("figures/Fig2.SplineCors.NN.pdf", width = 5, height = 9)

### MAP PROPORTIONS TO PRODUCE FIGURE 3

cities <- data.frame(Name = c("Cairns", "Townsville", "Gladstone"),
                     x = c(145.78, 146.82, 151.25),
                     y = c(-16.92, -19.26, -23.84))

# ReefsAUS.df <- read.csv("data/ReefsAUS.df.csv", strip.white = T)
# Very large spatial polygon file not included

Iso.plot <- ggplot() +
  geom_polygon(data = qld.spdf, aes(x = long, y = lat, group = group),
               fill = "#E9E8DF", colour = "darkgrey", size = .01) +
  geom_point(data = cities, aes(x = x, y = y), size = .5, shape = 4) +
  geom_text(data = cities, aes(x = x, y = y, label = Name),
            hjust = 1, nudge_x = -.2, size = 2) +
  # geom_polygon(data = ReefsAUS.df, aes(x = long, y = lat, group = group),
  #              fill = "#00F7DA", colour = NA) +
  geom_point(data = NN.emp, aes(x = Lon, y = Lat, colour = PropN*100), size=.3) +
  scale_colour_viridis_c(option = "A", direction = -1, end = .9, begin = .2) +
  facet_wrap(~ Event, ncol = 3) +
  coord_equal(xlim = c(142.5,153), ylim = c(-24,-8)) +
  scale_x_continuous(breaks = c(145,150)) +
  labs(x = "Longitude", y = "Latitude",
       colour = "% not severely disturbed\nreefs within 100 km radius") +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.background = element_rect(fill = "#B4C5CD"))

ggsave("figures/Fig3.MapPropN.pdf", width = 7, height = 6)
ggsave("figures/Fig3.MapPropN.png", width = 7, height = 6, dpi = 500)

### MAKE PAIRWISE COMPARISONS TO PRODUCE FIGURE S5

(Contrast.PropN.plot <- model.PropN.emp %>%
    emmeans(~ Event) %>%
    gather_emmeans_draws() %>%
    ungroup() %>%
    mutate(.value = inv_logit_scaled(.value)) %>%
    compare_levels(.value, by = Event) %>%
    ungroup() %>%
    mutate(Event = reorder(Event, .value)) %>%
    ggplot(aes(y = Event, x = .value)) +
    geom_halfeyeh(relative_scale = 2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Effect size", y = ""))

(Contrast.MinDist.plot <- model.MinDist.emp %>%
  emmeans(~ Event) %>%
  gather_emmeans_draws() %>%
  ungroup() %>%
  mutate(.value = (.value)) %>%
  compare_levels(.value, by = Event) %>%
  ungroup() %>%
  mutate(Event = reorder(Event, .value)) %>%
  ggplot(aes(y = Event, x = exp(.value))) +
  geom_halfeyeh() +
  scale_x_continuous(trans = "log", breaks = c(.5,1,2,4),
                     limits = c(.4,7)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Effect size", y = ""))

cowplot::plot_grid(Contrast.PropN.plot, Contrast.MinDist.plot, nrow = 2,
                   labels = c("(a)","(b)"), label_fontface = "italic")

ggsave("figures/FigS5.ContrastPlots.pdf", width = 6, height = 8)

### COMPARE ESTIMATES TO NULL EXPECTATIONS TO PRODUCE FIGURE S6

RelToNull <- as.data.frame(fits.NN[,c(1,6)]) %>%
  mutate(PropN = fits.PropN$PropN) %>%
  left_join(fits.MinDist.Null,
            by = "Event") %>%
  left_join(fits.PropN.Null,
            by = "Event") %>%
  mutate(NN.rel = exp(NN)/MinDist.Null,
         PropN.rel = PropN/PropN.Null)

RTN.PropN <- ggplot(RelToNull, aes(y = Event, x = 1/PropN.rel)) +
  geom_halfeyeh() +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_continuous(trans = "log", breaks = c(1,2,4,8,16,32)) +
  labs(x = "Factor\n(Proportion < null)", y = "")

RTN.MinDist <- ggplot(RelToNull, aes(y = Event, x = NN.rel)) +
  geom_halfeyeh() +
  geom_vline(xintercept = 1, linetype = 2) +
  theme(axis.text.y = element_blank()) +
  labs(x = "Factor\n(Distance > null)", y = "")

cowplot::plot_grid(RTN.PropN, RTN.MinDist, nrow = 1, rel_widths = c(1,.7),
                   labels = c("(a)","(b)"), label_fontface = "italic")

ggsave("figures/FigS6.RelativeToNull.pdf", height = 2.5, width = 7)

RelToNull %>% group_by(Event) %>%
  summarise(NN = median(NN.rel), PropN = 1/median(PropN.rel))

# SUMMARIZE FITS
fits.NN %>% group_by(Event) %>% summarise(NN = mean(exp(NN)))

fits.PropN %>% group_by(Event) %>% summarise(PropN = median(PropN))

#########################################################
# SENSITIVITY TO SEVERITY TRESHOLD ######################

# CALCULATE DISTANCE METRICS FOR EACH EVENT (SEVERITY THRESHOLD > 3)

c <- lapply(seq_along(unique(Events.Agg.df$Event)), function(i) {

  data = Events.Agg.df[Events.Agg.df$Event == unique(Events.Agg.df$Event)[i],] %>%
    mutate(ID1 = as.character(row_number()))

  coords = as.matrix(data.frame(Lon = data$Lon, Lat = data$Lat))

  distance.matrix <- as.data.frame(rdist.earth(coords, coords,
                                               miles = F, R = NULL))

  dist.mat.long <- distance.matrix %>%
    mutate(ID1 = as.character(1 : length(data$Lon))) %>%
    tidyr::gather(., ID2, distance, c(1 : length(data$Severe45))) %>%
    mutate(ID2 = substr(ID2, 2, 5),
           Self = as.numeric(ID1) - as.numeric(ID2)) %>%
    filter(!Self == 0) %>%
    left_join(., data[, c("ID1", "Severe45")], by = "ID1") %>%
    rename(Severe45A = Severe45) %>%
    left_join(., data[, c("ID1", "Severe45")], by = c("ID2" = "ID1")) %>%
    rename(Severe45B = Severe45)

  NN <- dist.mat.long %>%
    filter(Severe45B < 1) %>%
    filter(Severe45A > 0) %>%
    group_by(ID1) %>%
    summarise(MinDist = min(distance), .groups = "keep")

  PropN <- dist.mat.long %>% filter(distance < 100) %>%
    filter(Severe45A > 0) %>%
    group_by(ID1, Severe45B) %>%
    summarise(Freq = n(), .groups = "keep") %>% ungroup() %>%
    complete(ID1, Severe45B, fill = list(Freq = 0)) %>%
    group_by(ID1) %>%
    mutate(Prop = Freq/sum(Freq),
           FreqAll = sum(Freq)) %>%
    filter(Severe45B < 1) %>%
    dplyr::select(ID1, Prop, FreqAll, Freq) %>%
    left_join(., data[,c("Lon","Lat","ID1")], by="ID1")

  print(i)

  output = data.frame(Event = unique(Events.Agg.df$Event)[i],
                      MinDist = NN$MinDist,
                      PropN = PropN$Prop,
                      FreqAll = PropN$FreqAll,
                      Freq = PropN$Freq,
                      Lon = PropN$Lon,
                      Lat = PropN$Lat)
})

NN.emp.45 <- bind_rows(c) %>%
  mutate(type = "Empirical",
         severe = ">3")

### CREATE NULL MODEL DISTRIBUTIONS
### Shuffle column Binary n times
### Record for each event the mean distance and proportion

sims = 100

d <- lapply(seq_along(unique(Events.Agg.df$Event)), function(i) {

  data = Events.Agg.df[Events.Agg.df$Event == unique(Events.Agg.df$Event)[i], ] %>%
    mutate(ID1 = as.character(row_number()))

  coords = as.matrix(data.frame(Lon = data$Lon, Lat = data$Lat))

  distance.matrix <- as.data.frame(rdist.earth(coords, coords,
                                               miles = F, R = NULL))

  dist.mat.long <- distance.matrix %>%
    mutate(ID1 = as.character(1 : length(data$Lon))) %>%
    tidyr::gather(., ID2, distance, c(1 : length(data$Severe45))) %>%
    mutate(ID2 = substr(ID2, 2, 5),
           Self = as.numeric(ID1) - as.numeric(ID2)) %>%
    filter(!Self == 0
           # & distance > 0
    ) %>%
    #### Remove distance > 0 because Yasi scores have same coordinates
    left_join(., data[, c("ID1", "Severe45")], by = "ID1") %>%
    rename(Severe45A = Severe45) %>%
    left_join(., data[, c("ID1", "Severe45")], by = c("ID2" = "ID1")) %>%
    rename(Severe45B = Severe45)

  MinDist = rep(NA, sims)
  PropN = rep(NA, sims)

  for (j in 1 : sims) {
    dist.mat.long$Severe45B <- sample(dist.mat.long$Severe45B, replace = FALSE)

    NN <- dist.mat.long %>%
      filter(Severe45B < 1) %>%
      filter(Severe45A > 0) %>%
      group_by(ID1) %>%
      summarise(MinDist = min(distance), .groups = "keep")

    MinDist <- c(MinDist, NN$MinDist)

    PropNN <- dist.mat.long %>% filter(distance < 100) %>%
      filter(Severe45A > 0) %>%
      group_by(ID1, Severe45B) %>%
      summarise(Freq = n(), .groups = "keep") %>%
      ungroup() %>%
      complete(ID1, Severe45B, fill = list(Freq = 0)) %>%
      group_by(ID1) %>%
      mutate(Prop = Freq / sum(Freq)) %>%
      filter(Severe45B < 1) %>%
      dplyr::select(ID1, Prop) %>%
      left_join(., data[,c("Lon", "Lat", "ID1")], by = "ID1")

    PropN <- c(PropN, PropNN$Prop)
  }

  print(i)

  Null = data.frame(Event = unique(Events.Agg.df$Event)[i],
                    MinDist = MinDist,
                    PropN = PropN)
})

NN.Null.45 <- bind_rows(d) %>%
  filter(!is.na(MinDist)) %>%
  mutate(type = "Null", Lon = NA, Lat = NA,
         severe = ">3")


### MODEL FITS FOR NULL EXPECTATIONS

data.MinDist.45 <- NN.Null.45 %>%
  mutate(MinDist = ifelse(MinDist == 0, 0.01, MinDist))

model.MinDist.Null.45 <- brm(log(MinDist) ~ Event, data = data.MinDist.45)

fits.MinDist.Null.45 = NN.Null.45 %>%
  data_grid(Event) %>%
  add_fitted_draws(model.MinDist.Null.45) %>%
  mutate(MinDist = exp(.value)) %>%
  group_by(Event) %>%
  summarise(MinDist.Null = median(MinDist), .groups = "keep") %>%
  mutate(severe = ">3")

data.PropN.45 <- NN.Null.45 %>%
  mutate(PropN = ifelse(PropN == 0, PropN + 0.001,
                        ifelse(PropN == 1, PropN - 0.001, PropN)))

model.PropN.Null.45 <- brm(logit_scaled(PropN) ~ Event, data = data.PropN.45)

fits.PropN.Null.45 = NN.Null.45 %>%
  data_grid(Event) %>%
  add_fitted_draws(model.PropN.Null.45) %>%
  mutate(PropN = inv_logit_scaled(.value)) %>%
  group_by(Event) %>%
  summarise(PropN.Null = median(PropN), .groups = "keep") %>%
  mutate(severe = ">3")

### MODELS FITS FOR EMPIRICAL DATA

NN.emp2.45 <- NN.emp.45 %>%
  mutate(PropN = ifelse(PropN == 0, PropN + 0.001,
                        ifelse(PropN == 1, PropN - 0.001, PropN)))

model.PropN.emp.45 <- brm(bf(logit_scaled(PropN) ~ Event,
                          sigma ~ Event),
                       data = NN.emp2.45, family = "gaussian")

fits.PropN.45 <- NN.emp2.45 %>%
  data_grid(Event) %>%
  add_fitted_draws(model.PropN.emp.45) %>%
  rename(PropN = .value) %>%
  mutate(PropN = inv_logit_scaled(PropN),
         severe = ">3")

NN.emp2.45 <- NN.emp.45 %>%
  mutate(MinDist = ifelse(MinDist == 0, 0.01, MinDist))

model.MinDist.emp.45 <- brm(log(MinDist) ~ Event, data = NN.emp2.45,
                         family = "gaussian")

fits.NN.45 <- NN.emp2.45 %>%
  data_grid(Event) %>%
  add_fitted_draws(model.MinDist.emp.45) %>%
  rename(NN = .value) %>%
  mutate(severe = ">3")


### PLOTTING PROPORTIONS AND MINIMUM DISTANCES TO PRODUCE FIGURE 2

(ModelFit.PropN.plot <- bind_rows(NN.emp, NN.emp.45) %>%
    ggplot(aes(y = PropN, x = Event)) +
    geom_boxplot(aes(fill = severe), outlier.shape = NA) +
    geom_point(data = bind_rows(fits.PropN.Null, fits.PropN.Null.45),
               aes(x = Event, y = PropN.Null, group = severe),
               colour = "red", fill = "white", shape = 21, size = 2,
               position = position_dodge(width = .75)) +
    stat_pointinterval(data = bind_rows(fits.PropN, fits.PropN.45),
                       aes(y = PropN, group = severe),
                       .width = c(.66, .95),
                       position = position_dodge(width = .75), size = .5,
                       colour = "red") +
    labs(x = "", y = "Proportion", colour = "Interval") +
    scale_fill_manual(values = c("darkgrey", "white")) +
    guides(fill = F) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5),
          text = element_text(size = 18)))

(ModelFit.MinDist.plot <- bind_rows(NN.emp, NN.emp.45) %>%
    ggplot(aes(y = MinDist, x = Event)) +
    geom_boxplot(aes(fill = severe), outlier.shape = NA) +
    geom_point(data = bind_rows(fits.MinDist.Null, fits.MinDist.Null.45),
               aes(x = Event, y = MinDist.Null, group = severe),
               colour = "red", size = 2, shape = 21, fill = "white",
               position = position_dodge(width = .75)) +
    stat_pointinterval(aes(y = exp(NN), group = severe), colour = "red",
                       data = bind_rows(fits.NN, fits.NN.45),
                       .width = c(.66, .95),
                       position = position_dodge(width = .75), size = .5) +
    scale_y_continuous(limits = c(0,350)) +
    scale_fill_manual(values = c("darkgrey", "white")) +
    labs(x = "", y = "Distance (km)", colour = "Prediction\ninterval",
         fill = "Severity\nthreshold") +
    theme(axis.text.x = element_text(angle = 90, vjust = .5),
          legend.position = c(.25, .7),
          text = element_text(size = 18)))

SC.SensitivityThreshold <- SC.AggReef %>%
  group_by(event) %>%
  filter(distance < 0.5 * max(distance)) %>%
  ggplot(aes(x = distance, y = correlation,
             colour = severe, fill = severe)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(x = distance, ymin = lower, ymax = upper),
              alpha = .3, colour = NA) +
  geom_line() +
  facet_wrap(~ event, ncol = 2) +
  coord_cartesian(xlim = c(0, 350), ylim = c(-.3, 1)) +
  scale_colour_manual(values = c("black", "white")) +
  scale_fill_manual(values = c("grey", "black")) +
  labs(x = "Distance (km)", y = "Correlation",
       fill = "Severity threshold", colour = "Severity threshold") +
  theme(legend.position = "top", text = element_text(size = 18))


cowplot::plot_grid(SC.SensitivityThreshold,
                   cowplot::plot_grid(ModelFit.PropN.plot, ModelFit.MinDist.plot, nrow = 1,
                                      labels = c("(b)","(c)")),
                   ncol = 1, rel_heights = c(1.8, 1), labels = c("(a)",""))

ggsave("figures/FigSXX.SC_DistMetr_Sensitivity.pdf", height = 12, width = 7)

# SAVE OUTPUT FILES
save(NN.Null, NN.emp, model.PropN.emp, model.PropN.Null,
     model.MinDist.emp, model.MinDist.Null,
     fits.MinDist.Null, fits.PropN.Null,
     NN.Null.45, NN.emp.45, model.PropN.emp.45, model.PropN.Null.45,
     model.MinDist.emp.45, model.MinDist.Null.45,
     fits.MinDist.Null.45, fits.PropN.Null.45,
     fits.NN, fits.PropN, fits.NN.45, fits.PropN.45,
     RelToNull, file = "RData/3-NN.RData")
