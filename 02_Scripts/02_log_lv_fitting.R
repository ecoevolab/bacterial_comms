# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Logistic Lotka-Volterra fitting 

library(rstan)
library(readr)
library(dplyr)

# Reference --------------------------------------------------------
# Stan & Lotka-Volterra 
# https://canada1.discourse-cdn.com/flex030/uploads/mc_stan/original/2X/b/b6acad3d9f2cd8e7b925c90809f9fb4e4b005a4f.pdf 
 

# Data.frame --------------------------------------------------------
indv_gr <- read_tsv("01_RawData/individual_strains_growth_curves_filtered.tsv")
indv_gr <- as.data.frame(indv_gr)

# data.frame with aditional column for replica
indv_gro <- read_tsv("01_RawData/modified_individual_strain_growth_curves_ord_replica - individual_strains_growth_curves_filtered.tsv")
indv_gro <- as.data.frame(indv_gro)

# Pipeline -------------------------------------------------------------------
# we're going to use the data.frame without CH29 because the sampling is different from the other spps 

pruebafunc <- stan_ccfunct(df = indv_gro, temp_col = "temp", replica_col = "ord_replica", 
                           strain_col = "Cepa", sample_byh = "hr", interest_col = "OD_real", 
                           niterations = 1500, nchains = 3, rvalin = 0.3, kvalin = 0.8, sigma_val = 0.05)

rk_valslist_wtCH29 <- saveRDS(pruebafunc, file = "03_Output/stan_func_rk_vectors")
