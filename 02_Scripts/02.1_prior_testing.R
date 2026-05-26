# Author: 
# Title: Stan priors exploration 

library(tidyverse)
library(deSolve)

# Data -----------------------------------------------------------------------------------------
indv_gro <- read_tsv("01_RawData/modified_individual_strain_growth_curves_ord_replica - individual_strains_growth_curves_filtered.tsv")
indv_gro <- as.data.frame(indv_gro)

rslopedf <- read.csv("01_RawData/rk_slope_n_priors.csv")
rslopedf

# rk_slope_priors <- saveRDS(rslopedf, "03_Output/rk_slope_priors")
rk_slope_priors <- readRDS("03_Output/rkpriors_slope")


# r & k prior testing with heatmap results ---------------------------------------------------------------

rk_valz_cc <- rk_prior_testing(df = indv_gro, temp_col = "temp", replica_col = "ord_replica", strain_col = "Cepa", 
                               sample_byh = "hr", interest_col = "OD_real")
rk_valz_cc$grid_output





# r & k prior delimitation by data -------------------------------------------------------------
rk_slope_priors <- r_slope(df = indv_gro, xf = rslopedf, temp_col = "temp", strain_col = "Cepa", 
                    x1_col = "xin", x2_col = "xfin",samplebyh = "hr", interest_col = "OD_real")   
# new data.frame with slope column / DF YOU SHOULD USE TO TEST THE INITIAL VALUES 
