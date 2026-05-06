# Author: Mayra Beatriz Mendoza Velazquez
# Title: Cuatro ciénegas growth curves

library(readr)
library(ggplot2)
library(gridExtra)

# Load data.frame 
individual_growth_curves <- read_tsv("bacterial_comms/01_RawData/individual_strains_growth_curves_filtered.tsv")
individual_growth <- as.data.frame(individual_growth_curves)

#### Growth curves ####
gcurves_cc <- growth_curves_func(df = individual_growth, temps = c(30,37,42), strain_col = "Cepa", 
                                 temp_col = "temp", samplebyh = "hr", interest_column = "OD_real", 
                                 replica_col = "rep")

print(gcurves_cc)
