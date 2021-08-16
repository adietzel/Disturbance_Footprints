Data and code for Dietzel et al. 2021: The spatial footprint and patchiness of large-scale disturbances on coral reefs

Corresponding author: Andreas Dietzel, andreas.dietzel@my.jcu.edu.au

Abstract:

Ecosystems have always been shaped by disturbances, but many of these events are becoming larger, more severe and more frequent. The recovery capacity of depleted populations depends on the frequency of disturbances, the spatial distribution of mortality and the scale of dispersal. Here, we show that four mass coral bleaching events on the Great Barrier Reef (in 1998, 2002, 2016 and 2017) each had markedly larger disturbance footprints and were less patchy than a severe category 5 tropical cyclone (Cyclone Yasi, 2011). Severely bleached reefs in 2016 and 2017 were isolated from the nearest lightly affected reefs by up to 146 and 200 km, respectively. In contrast, reefs damaged by Cyclone Yasi were on average 20 km away from relatively undisturbed reefs, well within the estimated range of larval dispersal for most corals. Based on these results, we present a model of coral reef disturbance and recovery to examine (1) how the spatial clustering of disturbances modifies large-scale recovery rates; and (2) how recovery rates are shaped by species' dispersal abilities. Our findings illustrate that the spatial footprint of the recent mass bleaching events poses an unprecedented threat to the resilience of coral species in human history, a threat that is even larger than the amount of mortality suggests.

Folder structure:

The file "Spatial_Footprint_Disturbances.Rproj" sets the directory and work environment. 
R scripts can be found in subfolder "R", numbered in sequence. 
Input data can be found in subfolder "data". 
Figures are automatically stored in the subfolder "figures".
Files saving intermediate results are stored in the subfolder "RData".
MatLab scripts to run reef disturbance and recovery model are stored in the subfolder "Matlab"

Input files:

File "data/Disturbance_Events.csv"
Data of impact of individual disturbance events.

File "RData/0-Auxiliary.RData"
Important auxiliary data such as reef area, map of Queensland for plotting, spatial grid.

File "data/TropCyclones.xlsx" 
Data showing tracks of cyclones crossing the Great Barrier Reef.

File "data/Cyclones_data_updated2017.csv"
Additional data showing disturbance histories of individual reefs on the Great Barrier Reef. Used to calculate cumulative hours of destructive wind speeds (Figure 1). From Matthews et al. 2019 (https://doi.org/10.1002/ecy.2574). 

File "data/Centroids.csv"
Locations of centroids of individual reefs on the Great Barrier Reef. Used to generate spatial disturbance patterns. 
