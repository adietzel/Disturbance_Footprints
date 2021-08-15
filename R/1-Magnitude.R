### The magnitude (extent x severity) of reef disturbance events ######
### Author: Andreas Dietzel, andreas.dietzel@my.jcu.edu.au ############


rm(list = ls())

# LOAD LIBRARIES AND DATA
library(tidyverse)
library(cowplot)
library(raster)
library(sp)

load("RData/0-Auxiliary.RData")
load("RData/0-Events.RData")
theme_set(theme_simple)

# CALCULATE REEF AREA IN EACH GRID CELL

Rast.events.df <- data.frame(
  Cell = 1 : length(reefs.rast.agg@data@values),
  # Proportion of reef area in each grid cell times size in km2
  ReefArea.km2 = reefs.rast.agg@data@values * 50 ^ 2)

# IDENTIFY IN WHICH GRID CELL EACH OBSERVATION IS LOCATED

Events.df$Cell <- over(Events.cea, as(reefs.rast.agg, "SpatialPixels"))

# CALCULATE TOTAL REEF AREA ON THE GBR
TotalReefArea.GBR <- sum(reefs.rast.agg@data@values *
                           50 ^ 2, na.rm = T)

# CALCULATE CUMULATIVE % BLEACHED CORALS IN EACH GRID CELL ACROSS THE BACK-TO-BACK
# BLEACHING EVENT AS
# AVERAGE % BLEACHED IN 2016 + % NOT BLEACHED IN 2016 * % BLEACHED IN 2017

Magnitude.1617 <- Events.df %>%
  filter(Event %in% c("2016 Bleaching", "2017 Bleaching")) %>%
  mutate(Event = ifelse(Event == "2016 Bleaching", "Bleach16", "Bleach17")) %>%
  group_by(Event, Cell) %>%
  summarise(IntMean = mean(IntMean), .groups = "keep") %>%
  mutate(IntMean = ifelse(is.na(IntMean), 0 , IntMean)) %>%
  spread(key = Event, value = IntMean) %>%
  replace_na(., list(Bleach16 = 0, Bleach17 = 0)) %>%
  mutate(IntMean = Bleach16 + ((100 - Bleach16) * Bleach17/100),
         Level = ifelse(
           IntMean < 5, 0,
           ifelse(IntMean >= 5 & IntMean < 10, 1,
                  ifelse(IntMean >= 10 & IntMean < 30, 2,
                         ifelse(IntMean >= 30 & IntMean < 60, 3, 4)))),
         Event = "Back-to-back") %>%
  dplyr::select(-Bleach16, -Bleach17)

# CACLULATE MAGNITUDE OF EACH EVENT AS % OF GBR REEF AREA AFFECTED BY DIFFERENT
# SEVERITIES

Magnitude.Events <- Events.df %>%
  # append back-to-back magnitude data calculate above
  bind_rows(., Magnitude.1617) %>%
  # join reef area data
  full_join(Rast.events.df, by = "Cell") %>%
  # exclude cells with no reef area
  filter(ReefArea.km2 > 0) %>%
  # ensure each disturbance has all grid cells
  complete(Event, nesting(Cell, ReefArea.km2)) %>%
  group_by(Event) %>%
  # if grid cell with no severity score, assign score of 0
  mutate(Level = ifelse(is.na(Level), 0, Level)) %>%
  group_by(Event, Cell) %>%
  # distribute reef area in each grid cell equally across observations
  mutate(ReefArea.prop = ReefArea.km2/length(ReefArea.km2)) %>%
  group_by(Event, Level) %>%
  # calculate for each event total reef area with differen severity
  summarise(TotReefArea = sum(ReefArea.prop), .groups = "keep") %>%
  ungroup() %>%
  # calculate % of GBR reef area
  mutate(PropGBR = TotReefArea/TotalReefArea.GBR,
         Event = as.factor(Event),
         Event = fct_reorder(Event, desc(Event))) %>%
  filter(Level > 0) %>%
  mutate(Level = ifelse(Level > 3, "4-5", Level))

# PLOT THE RESULTS

(Magnitude.plot <- ggplot(
  Magnitude.Events, aes(x = forcats::fct_rev(Event),
                        y = PropGBR*100, fill = as.factor(Level))) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d(option = "D") +
    labs(y = "Percentage of GBR", x = "", fill = "Severity") +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 90, vjust = .5)) +
    ylim(0, 100))

ggsave("figures/Fig1b.MagnitudePlot.pdf", width = 2.5, height = 3.5)

# BRIEF SUMMARY TABLE OF RESULTS

Magnitude.Events %>%
  mutate(PropGBR = round(PropGBR * 100, digits = 1)) %>%
  dplyr::select(-TotReefArea) %>%
  pivot_wider(names_from = Event, values_from = PropGBR)

save(Magnitude.Events, file = "RData/1-Magnitude.RData")

# CUMULATIVE BLEACHING VS CYCLONES

CumulativeImpact <- Magnitude.Events %>%
  mutate(Type = ifelse(Event == "2011 Cyclone Yasi", "Cyclone Yasi", "Bleaching")) %>%
  group_by(Type, Level) %>%
  summarise(Cum_ReefArea = sum(TotReefArea),
            Cum_PropGBR = sum(PropGBR),
            .groups = "keep")

CumImpactYasi23 <- CumulativeImpact %>%
  filter(Type == "Cyclone Yasi") %>%
  mutate(Type = "Cyclone Yasi * 23",
         Cum_ReefArea = Cum_ReefArea * 23,
         Cum_PropGBR = Cum_PropGBR * 23)

bind_rows(CumulativeImpact, CumImpactYasi23) %>%
  ggplot(aes(x = Type, y = Cum_PropGBR, fill = as.factor(Level))) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d(option = "D") +
    labs(y = "Proportion of GBR", x = "", fill = "Severity") +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))

ggsave("figures/CumulativeImpactComparison.pdf", width = 2.3, height = 3)


## FIGURE 1

cities <- data.frame(Name = c("Cairns", "Townsville", "Gladstone"),
                     x = c(145.78, 146.82, 151.25),
                     y = c(-16.92, -19.26, -23.84))

# ReefsAUS.df <- read.csv("data/ReefsAUS.df.csv", strip.white = T)
# Very large spatial polygon file not included

Events.df2 <- Events.df %>%
  filter(Event %in% c("1998 Bleaching Agg Reef", "2002 Bleaching Agg Reef",
                      "2011 Yasi Agg Reef", "2016 Bleaching Agg Reef",
                      "2017 Bleaching Agg Reef")) %>%
  mutate(Event = gsub(" Agg Reef", "", Event)) %>%
  mutate(Event = gsub("Yasi", "Cyclone Yasi", Event))

FootprintsBinary.plot <- ggplot() +
  geom_polygon(data = qld.spdf, aes(x = long, y = lat, group = group),
               fill = "#E9E8DF", colour = "darkgrey", size = .01) +
  geom_point(data = cities, aes(x = x, y = y), size = .5, shape = 4) +
  geom_text(data = cities, aes(x = x, y = y, label = Name),
            hjust = 1, nudge_x = -.2, size = 2) +
  # geom_polygon(data = ReefsAUS.df, aes(x = long, y = lat, group = group),
  #              fill = "#00F7DA", colour = NA) +
  geom_point(data = Events.Agg.df[!Events.Agg.df$Event == "Back-to-back",],
             aes(x = Lon, y = Lat, colour = as.factor(Severe345)),
             size = .5, shape = 16, alpha = .7) +
  scale_colour_viridis_d(option = "A", end = .8,
                         labels = c("not severe","severe")) +
  facet_wrap(~ Event, ncol = 3) +
  coord_equal(xlim = c(142.5,153), ylim = c(-24,-8)) +
  scale_x_continuous(breaks = c(145,150)) +
  labs(x = "Longitude", y = "Latitude",
       colour = "Severity") +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.background = element_rect(fill = "#B4C5CD"))

ggsave("figures/Fig1a.Footprints.pdf", width = 8, height = 6)

### Tracks of cyclones since 1998
library(readxl)
Cyclones <- read_excel("data/TropCyclones.xlsx", trim_ws = T,
                       col_types = c(rep("guess",4),rep("numeric",3)))

Cyclones.sel <- Cyclones %>%
  mutate(Year = as.numeric(substr(TM, 0, 4))) %>%
  filter(Year > 1997) %>%
  mutate(Category = ifelse(
    MAX_WIND_SPD < 118/3.6, "< Cat 3",
    ifelse(MAX_WIND_SPD > 118/3.6 & MAX_WIND_SPD < 158/3.6, "Cat 3",
           ifelse(MAX_WIND_SPD > 158/3.6 & MAX_WIND_SPD < 198/3.6, "Cat 4",
                  "Cat 5")))
  ) %>%
  arrange(TM) %>%
  # filter(!(is.na(CycloneT) | CycloneT == "< Cat3")) %>%
  filter(LON > 142.5 & LON < 152.5) %>%
  filter(LAT > -23.5 & LAT < -9.8)

Cyclones.sel.summary <- Cyclones.sel %>%
  group_by(NAME, Year) %>%
  summarise(MAX_WIND_SPD = max(MAX_WIND_SPD, na.rm = T),
            .groups = "keep") %>%
  mutate(Category = ifelse(
    MAX_WIND_SPD < 118/3.6, "< Cat 3",
    ifelse(MAX_WIND_SPD > 118/3.6 & MAX_WIND_SPD < 158/3.6, "Cat 3",
           ifelse(MAX_WIND_SPD > 158/3.6 & MAX_WIND_SPD < 198/3.6, "Cat 4",
                  "Cat 5")))) %>%
  filter(MAX_WIND_SPD > 62/3.6) %>%
  filter(!Category == "< Cat 3") %>%
  arrange(-Year) %>%
  filter(!NAME %in% c("Guba","KATRINA"))

Cyclones.df <- data.frame()

for (i in 1:length(unique(Cyclones.sel$NAME))) {
  plot.new()
  CycloneSel.ps <- data.frame(
    xspline(Cyclones.sel[Cyclones.sel$NAME == unique(Cyclones.sel$NAME)[i],5:6],
            shape=.8, lwd=2, draw=F),
    NAME = unique(Cyclones.sel$NAME)[i])
  Cyclones.df <- bind_rows(Cyclones.df, CycloneSel.ps)
}

CycloneData <- read.csv("data/Cyclones_data_updated2017.csv") %>%
  tidyr::gather(key, value, 7:39) %>%
  mutate(Year = substr(key, 7,11)) %>%
  group_by(REEF_ID, Year, LAT, LONG) %>%
  summarise(Damage = mean(value), .groups = "keep") %>%
  filter(Year > 1997) %>%
  group_by(REEF_ID, LAT, LONG) %>%
  summarise(Damage = sum(Damage), .groups = "keep")

CycloneHistory.plot <- ggplot() +
  geom_polygon(data = qld.spdf, aes(x = long, y = lat, group = group),
               fill = "#E9E8DF", colour = "darkgrey", size = .01) +
  geom_point(data = cities, aes(x = x, y = y), size = .5, shape = 4) +
  geom_text(data = cities, aes(x = x, y = y, label = Name),
            hjust = 1, nudge_x = -.2, size = 2) +
  geom_point(data = CycloneData[CycloneData$Damage > 0, ],
             aes(x = LONG, y = LAT, colour = Damage), size = .2) +
  geom_path(data = Cyclones.df[Cyclones.df$NAME %in% Cyclones.sel.summary$NAME,],
            colour = "black", aes(x = y, y = x, by = NAME)) +
  theme_simple +
  coord_equal(xlim = c(142.5,153), ylim = c(-24,-8)) +
  scale_x_continuous(breaks = c(145,150)) +
  labs(x = "Longitude", y = "Latitude", fill = "Cyclone category",
       colour = "Hours of\nhigh wind") +
  scale_colour_viridis_c(option = "A", end = .9, begin = .2) +
  scale_fill_viridis_d(begin = .3) +
  theme(axis.title = element_text(size = 14),
        panel.background = element_rect(fill = "#B4C5CD"))

# COMBINE PLOTS TO CREATE FIGURE 1

Figure <- cowplot::plot_grid(
  FootprintsBinary.plot, ncol = 1, rel_heights = c(1.5, 1), labels = c("(a)",""),
  cowplot::plot_grid(Magnitude.plot, CycloneHistory.plot,
                     nrow = 1, labels = c("(b)","(c)"), align = "h",
                     rel_widths = c(1,1.3)))

ggsave("figures/Fig1.FootprintsMagnitudeCycloneHistory.pdf", width = 6, height = 9)
