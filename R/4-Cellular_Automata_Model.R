### Analysing and plotting matlab output for conceptual model #########
### Author: Andreas Dietzel, andreas.dietzel@my.jcu.edu.au ############

rm(list = ls())

# LOAD PACKAGES AND DATA

library(tidyverse)
library(reshape2)

IB_High <- read.csv("Matlab/IB_1.csv", strip.white = T, header = F)
IB_Moderate <- read.csv("Matlab/IB_2.csv", strip.white = T, header = F)
IB_Low <- read.csv("Matlab/IB_3.csv", strip.white = T, header = F)
IB_No <- read.csv("Matlab/IB_4.csv", strip.white = T, header = F)

IS_High <- read.csv("Matlab/IS_1.csv", strip.white = T, header = F)
IS_Moderate <- read.csv("Matlab/IS_2.csv", strip.white = T, header = F)
IS_Low <- read.csv("Matlab/IS_3.csv", strip.white = T, header = F)
IS_No <- read.csv("Matlab/IS_4.csv", strip.white = T, header = F)

# LOAD AND SET CUSTOM THEME
load("RData/0-Events.RData")
load("RData/0-Auxiliary.RData")
theme_set(theme_simple)

### CALCULATE IMPACT OF SIMULATIONS

n = 10

Impact_Brood_High <- IB_High %>%
  mutate(Magnitude = seq(0, 1, length.out = 15),
         Mode = "Short-distance", SAC = "High") %>%
  tidyr::gather(key = rep, value = Impact, 1:n)

Impact_Brood_Moderate <- IB_Moderate %>%
  mutate(Magnitude = seq(0, 1, length.out = 15),
         Mode = "Short-distance", SAC = "Moderate") %>%
  tidyr::gather(key = rep, value = Impact, 1:n)

Impact_Brood_Low <- IB_Low %>%
  mutate(Magnitude = seq(0, 1, length.out = 15),
         Mode = "Short-distance", SAC = "Low") %>%
  tidyr::gather(key = rep, value = Impact, 1:n)

Impact_Brood_No <- IB_No %>%
  mutate(Magnitude = seq(0, 1, length.out = 15),
         Mode = "Short-distance", SAC = "None") %>%
  tidyr::gather(key = rep, value = Impact, 1:n)

Impact_Spawn_High <- IS_High %>%
  mutate(Magnitude = seq(0, 1, length.out = 15),
         Mode = "Long-distance", SAC = "High") %>%
  tidyr::gather(key = rep, value = Impact, 1:n)

Impact_Spawn_Moderate <- IS_Moderate %>%
  mutate(Magnitude = seq(0, 1, length.out = 15),
         Mode = "Long-distance", SAC = "Moderate") %>%
  tidyr::gather(key = rep, value = Impact, 1:n)

Impact_Spawn_Low <- IS_Low %>%
  mutate(Magnitude = seq(0, 1, length.out = 15),
         Mode = "Long-distance", SAC = "Low") %>%
  tidyr::gather(key = rep, value = Impact, 1:n)

Impact_Spawn_No <- IS_No %>%
  mutate(Magnitude = seq(0, 1, length.out = 15),
         Mode = "Long-distance", SAC = "None") %>%
  tidyr::gather(key = rep, value = Impact, 1:n)

Impact_all <- bind_rows(Impact_Brood_High, Impact_Brood_Moderate,
                        Impact_Brood_Low, Impact_Brood_No,
                        Impact_Spawn_High, Impact_Spawn_Moderate,
                        Impact_Spawn_Low, Impact_Spawn_No)

Impact_all$SAC <- as.factor(Impact_all$SAC)
Impact_all$SAC <- factor(Impact_all$SAC, levels(Impact_all$SAC)[c(1,3,2,4)])

# CALCULATE IMPACT AMPLIFICATION POTENTIAL

Impact_Amp_B <- Impact_all %>%
  filter(Mode == "Short-distance") %>%
  spread(key = SAC, value = Impact) %>%
  mutate(AmpHN = High/None, AmpMN = Moderate/None, AmpLN = Low/None) %>%
  dplyr::select(Mode, Magnitude, AmpHN, AmpMN, AmpLN) %>%
  tidyr::gather(key = Comp, value = Amp, 3:5)

Impact_Amp_S <- Impact_all %>%
  filter(Mode == "Long-distance") %>%
  spread(key = SAC, value = Impact) %>%
  mutate(AmpHN = High/None, AmpMN = Moderate/None, AmpLN = Low/None) %>%
  dplyr::select(Mode, Magnitude, AmpHN, AmpMN, AmpLN) %>%
  tidyr::gather(key = Comp, value = Amp, 3:5)

Impact_Amp <- bind_rows(Impact_Amp_B, Impact_Amp_S) %>%
  mutate(Comp = as.factor(Comp))

levels(Impact_Amp$Comp) <- c("High:None", "Low:None", "Moderate:None")

Impact_Amp$Comp <- factor(Impact_Amp$Comp, levels(Impact_Amp$Comp)[c(1,3,2)])

### PLOTTING RESULTS

(SAC.comp.plot <- ggplot(Impact_all, aes(x=Magnitude, y=Impact, colour=SAC)) +
    geom_smooth() +
    facet_wrap(~Mode) +
    scale_x_continuous(breaks=c(0,.5,1),labels=c("0","0.5","1")) +
    scale_colour_viridis_d(option="D") +
    theme(axis.line = element_blank()) +
    labs(colour="Spatial\nautocorrelation",
         x = "Magnitude (proportion of cells disturbed)",
         y = expression(paste("Impact ( ", integral(Delta),"N"[t], " dt )"))))

(Impact_Amp_Mode.plot <- Impact_all %>%
  spread(key = Mode, value = Impact) %>%
  mutate(Ratio = `Short-distance`/`Long-distance`) %>%
  ggplot(aes(x = Magnitude, y = Ratio, colour = SAC)) +
  theme(axis.line = element_blank()) +
  geom_smooth() +
  scale_x_continuous(breaks = c(0,.5,1), labels = c("0","0.5","1")) +
  scale_colour_viridis_d(option = "D") +
  labs(colour="Spatial\nautocorrelation",
       x = "Magnitude (proportion of cells disturbed)",
       y = expression(paste("Impact"[SD], " / Impact"[LD]))))

cowplot::plot_grid(
  cowplot::plot_grid(SAC.comp.plot + theme(legend.position = ""),
                     Impact_Amp_Mode.plot + ggtitle(""), align = "v",
                     ncol = 2, labels = c("(a)","(b)"), rel_widths = c(1,.8),
                     label_fontface = "bold"))

# ggsave("figures/Fig5ab.ModelOutput.pdf", width = 9, height = 2.5)

(Amplif.plot <-
    ggplot(Impact_Amp, aes(x = Magnitude, y = Amp, colour = Comp)) +
    geom_hline(yintercept = 1, linetype = 2) +
    theme(axis.line = element_blank()) +
    facet_wrap(~ Mode) +
    labs(y = "Amplification", colour = "Comparison") +
    geom_smooth() +
    labs(x = "Magnitude (proportion of cells disturbed)",
         y = expression(paste("Impact"["H,M,L"], " / Impact"[N]))) +
    scale_x_continuous(breaks = c(0,.5,1), labels = c("0","0.5","1")) +
    coord_cartesian(ylim = c(.95, 2.1)))

# ggsave("figures/Fig5c.ModelOutput.pdf", width = 6.7, height = 2.5)

# PLOT DISPERSAL KERNELS

Spawner <- matrix(data = c(rep(0.055/16, 5),
                           0.055/16, rep(0.035/8,3), 0.055/16,
                           0.055/16, 0.035/8, 0.01, 0.035/8, 0.055/16,
                           0.055/16, rep(0.035/8,3), 0.055/16,
                           rep(0.055/16, 5)),
                  nrow = 5, ncol = 5)

colnames(Spawner) <- 1:5
rownames(Spawner) <- 1:5

Brooder <- matrix(data = c(rep(NA, 5),
                           NA, rep(0.05/8,3), NA,
                           NA, 0.05/8, 0.05, 0.05/8, NA,
                           NA, rep(0.05/8,3), NA,
                           rep(NA, 5)),
                  nrow = 5, ncol = 5)

colnames(Brooder) <- 1:5
rownames(Brooder) <- 1:5

Spawn.long <- melt(Spawner) %>% mutate(type = "Long-distance")
Brood.long <- melt(Brooder) %>% mutate(type = "Short-distance")

Long <- bind_rows(Spawn.long, Brood.long)

KernelMat <- ggplot(Long, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value/.2*100)) +
  facet_wrap(~ type, ncol = 2) +
  scale_fill_viridis_c(option = "D") +
  labs(x = "column", y = "row", fill = "% larvae") +
  scale_x_continuous(breaks = c(1:5), labels = c(-2,-1,0,1,2)) +
  scale_y_continuous(breaks = c(1:5), labels = c(-2,-1,0,1,2)) +
  coord_equal()

KernelBar <- data.frame(Type = c(rep("Long-distance",3),rep("Short-distance",3)),
           Neighbour = rep(c("Local retention","Adjacent","Second order"),2),
           Value = c(.01,.035,.055,.05,.05,0)/0.1,
           Order = c(1:3, 1:3)) %>%
  ggplot(aes(x = reorder(Neighbour, Order), y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Proportion of larvae", fill = "")

cowplot::plot_grid(KernelMat, KernelBar, ncol = 1, align = "v",
                   rel_heights = c(1,.7))

ggsave("figures/FigS8.DispersalKernels.pdf", width = 5, height = 6)
