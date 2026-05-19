# Author: 
# Title: Stan priors exploration 

library(tidyverse)
library(deSolve)

# Data -------------------------------------------------------------------------
indv_gro <- read_tsv("01_RawData/modified_individual_strain_growth_curves_ord_replica - individual_strains_growth_curves_filtered.tsv")
indv_gro <- as.data.frame(indv_gro)

# r & k prior testing ---------------------------------------------------------------

rk_valz_cc <- rk_prior_testing(df = indv_gro, temp_col = "temp", replica_col = "ord_replica", strain_col = "Cepa", 
                               sample_byh = "hr", interest_col = "OD_real")

rk_valz_cc$grid_output$CH161d_T30

# saveRDS(rk_valz_cc, file = "03_Output/rkpriors_heatmap")
# rk_valz_cc<-  readRDS("03_Output/rkpriors_heatmap")
