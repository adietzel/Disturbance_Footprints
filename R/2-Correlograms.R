### PRODUCE CORRELOGRAMS FOR EACH DISTURBANCE EVENT ###################
### Author: Andreas Dietzel, andreas.dietzel@my.jcu.edu.au ############

rm(list = ls())

# LOAD LIBRARIES AND DATA
library(tidyverse)
require(ncf)

load("RData/0-Events.RData")
load("RData/0-Auxiliary.RData")

theme_set(theme_simple)

### PRODUCE CORRELOGRAMS FOR REEF AGGREGATE DATA TO PRODUCE FIGURE 2

Events.Agg.df <- Events.Agg.df %>%
  filter(Core %in% c("1","DW","VDW"))

SC.AggReef = data.frame()

# Produce correlogram for each event in for loop and store output
for (i in 1:length(unique(Events.Agg.df$Event))) {

  Event = unique(Events.Agg.df$Event)[i]
  print(Event)
  data = Events.Agg.df[Events.Agg.df$Event == Event,]

  SC345 <- spline.correlog(data$Lon, data$Lat, data$Severe345,
                        npoints = 300, resamp = 100, latlon = T)

  df345 <- data.frame(event = Event,
                   distance = t(SC345$real$predicted$x),
                   correlation = t(SC345$real$predicted$y),
                   lower = t(SC345$boot$boot.summary$predicted$y)[,2],
                   upper = t(SC345$boot$boot.summary$predicted$y)[,10],
                   type = "Reef aggregate",
                   severe = ">2")

  SC45 <- spline.correlog(data$Lon, data$Lat, data$Severe45,
                        npoints = 300, resamp = 100, latlon = T)

  df45 <- data.frame(event = Event,
                   distance = t(SC45$real$predicted$x),
                   correlation = t(SC45$real$predicted$y),
                   lower = t(SC45$boot$boot.summary$predicted$y)[,2],
                   upper = t(SC45$boot$boot.summary$predicted$y)[,10],
                   type = "Reef aggregate",
                   severe = ">3")

  SC.AggReef = bind_rows(SC.AggReef, df345, df45)

}

SC.AggReef %>% group_by(event) %>%
  filter(distance < 0.5 * max(distance)) %>%
  # filter(severe == ">2") %>%
  ggplot(aes(x = distance, y = correlation,
             colour = severe, fill = severe)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(x = distance, ymin = lower, ymax = upper),
              alpha = .3, colour = NA) +
  geom_line() +
  facet_wrap(~ event, ncol = 2) +
  coord_cartesian(xlim = c(0, 350), ylim = c(-.3, 1)) +
  labs(x = "Distance (km)", y = "Correlation")

ggsave("figures/FigSXX.SplineCors.SensitivityThreshold.pdf", width = 4, height = 5)


#### SENSITIVITY OF CORRELOGRAMS TO SCALE TRANSFORMATION
#### Ordinal VS binary VS interval mid-points

Events.df.sel <- Events.df %>%
  filter(!is.na(Level)) %>%
  filter(Core %in% c("1","DW","VDW"))

SC.ScaleComp = data.frame()

for (i in 1:length(unique(Events.df.sel$Event))) {

  Event = unique(Events.df.sel$Event)[i]
  print(Event)

  data = Events.df.sel[Events.df.sel$Event == Event,]

  SC.Site <- spline.correlog(data$Lon, data$Lat, data$Binary,
                             npoints = 300, resamp = 100, latlon = T)
  SC.Ord <- spline.correlog(data$Lon, data$Lat, data$Level,
                            npoints = 300, resamp = 100, latlon = T)
  SC.IntMean <- spline.correlog(data$Lon, data$Lat, data$IntMean,
                                npoints = 300, resamp = 100, latlon = T)

  df.Site <- data.frame(event = Event,
                        distance = t(SC.Site$real$predicted$x),
                        correlation = t(SC.Site$real$predicted$y),
                        lower = t(SC.Site$boot$boot.summary$predicted$y)[,2],
                        upper = t(SC.Site$boot$boot.summary$predicted$y)[,10],
                        type = "Binary")

  df.Ord <- data.frame(event = Event,
                       distance = t(SC.Ord$real$predicted$x),
                       correlation = t(SC.Ord$real$predicted$y),
                       lower = t(SC.Ord$boot$boot.summary$predicted$y)[,2],
                       upper = t(SC.Ord$boot$boot.summary$predicted$y)[,10],
                       type = "Ordinal")

  df.IntMean <- data.frame(event = Event,
                           distance = t(SC.IntMean$real$predicted$x),
                           correlation = t(SC.IntMean$real$predicted$y),
                           lower = t(SC.IntMean$boot$boot.summary$predicted$y)[,2],
                           upper = t(SC.IntMean$boot$boot.summary$predicted$y)[,10],
                           type = "Interval mid-point")

  df.temp = bind_rows(df.Site, df.Ord, df.IntMean)

  SC.ScaleComp = bind_rows(SC.ScaleComp, df.temp)

}

SC.ScaleComp %>%
  group_by(event) %>%
  filter(distance < .5 * max(distance)) %>%
  ggplot(aes(x = distance, y = correlation, linetype = type)) +
  geom_line() +
  facet_wrap(~ event, ncol = 2) +
  coord_cartesian(xlim = c(0, 350), ylim = c(-.3, 1)) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(linetype = "", x = "Distance (km)", y = "Correlation")

ggsave("figures/FigS3.SplineCor.ScaleComp.pdf", width = 5, height = 5)

### STORE OUTPUT IN RDATA-FILE

save(SC.AggReef, SC.ScaleComp, file = "RData/2-SplineCor.RData")
