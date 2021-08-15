### Preprocessing data ################################################
### Author: Andreas Dietzel, andreas.dietzel@my.jcu.edu.au ############

rm(list = ls())

# Load libraries and files ####
library(tidyverse)

### LOAD DATA
Events.df <- read.csv("data/Disturbance_Events.csv", strip.white = T)
load("RData/0-Auxiliary.RData") # spatial grid, map of Queensland, raster with reef area ...

Events.df$Binary <- ifelse(Events.df$Level > 2, 1, 0)

### CALCULATE REEF AVERAGES

YasiAggReef <- Events.df %>%
  filter(Event == "2011 Cyclone Yasi") %>%
  group_by(Event, Lat, Lon, ReefID, Core) %>%
  summarise(Level_Matrix = mean(Level_Matrix), .groups = "keep") %>%
  mutate(Severe345 = ifelse(Level_Matrix > 90, 1, 0),
         Severe45 = ifelse(Level_Matrix > 135, 1, 0),
         Level = NA, IntMean = NA) %>%
  ungroup() %>% as.data.frame()

BleachAggReef <- Events.df %>%
  filter(!Event == "2011 Cyclone Yasi") %>%
  group_by(Event, ReefID, Core) %>%
  summarise(IntMean = mean(IntMean),
            Lat = mean(Lat),
            Lon = mean(Lon),
            .groups = "keep") %>%
  mutate(Severe345 = ifelse(IntMean > 30, 1, 0),
         Severe45 = ifelse(IntMean > 60, 1, 0),
         Level = NA, Level_Matrix = NA) %>%
  ungroup() %>% as.data.frame()

# BACK-TO-BACK BLEACHING EVENT

Reefs1617 <- BleachAggReef %>%
  filter(Event %in% c("2016 Bleaching", "2017 Bleaching")) %>%
  group_by(ReefID) %>%
  summarise(Freq = n(), .groups = "keep") %>%
  filter(Freq == 2)

BackToBack <- BleachAggReef %>%
  filter(Event %in% c("2016 Bleaching", "2017 Bleaching") &
           ReefID %in% Reefs1617$ReefID) %>%
  group_by(Event, ReefID) %>%
  summarise(IntMean = mean(IntMean),
            Lat = mean(Lat),
            Lon = mean(Lon),
            .groups = "keep") %>%
  mutate(IntMean = ifelse(is.na(IntMean), 0 , IntMean)) %>%
  spread(key = Event, value = IntMean) %>%
  replace_na(., list(`2016 Bleaching` = 0, `2017 Bleaching` = 0)) %>%
  group_by(ReefID) %>%
  summarise(Lat = mean(Lat), Lon = mean(Lon),
            Bleach16 = sum(`2016 Bleaching`),
            Bleach17 = sum(`2017 Bleaching`),
            .groups = "keep") %>%
  mutate(IntMean = Bleach16 + ((100 - Bleach16) * Bleach17/100),
         Level = ifelse(IntMean < 1, 0,
                        ifelse(IntMean >= 1 & IntMean < 10, 1,
                               ifelse(IntMean >= 10 & IntMean < 30, 2,
                                      ifelse(IntMean >= 30 & IntMean < 60, 3, 4))))) %>%
  dplyr::select(-Bleach16, -Bleach17) %>%
  mutate(Severe345 = ifelse(IntMean > 30, 1, 0),
         Severe45 = ifelse(IntMean > 60, 1, 0),
         Event = "Back-to-back",
         Level_Matrix = NA,
         Core = ifelse(Lat > -22, 1, 0)) %>%
  ungroup() %>% as.data.frame()

Events.Agg.df <- rbind(BleachAggReef, YasiAggReef, BackToBack)

# CONVERT TO SPATIAL POINTS DATA FRAME
Events.spdf <- SpatialPointsDataFrame(coords = Events.df[, c("Lon", "Lat")],
                                      data = Events.df,
                                      proj4string = CRS(proj_LatLong))

Events.cea <- spTransform(Events.spdf, CRS(proj_CEA))

# SAVE OUTPUT FILES

save(Events.df, Events.Agg.df, Events.cea, file = "RData/0-Events.RData")

