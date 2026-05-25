# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Logistic Lotka-Volterra fitting 

library(rstan)
library(readr)
library(dplyr)
library(deSolve)
library(tidytable)
library(fitdistrplus)
library(ggplot2)

# Reference --------------------------------------------------------
# Stan & Lotka-Volterra 
# https://canada1.discourse-cdn.com/flex030/uploads/mc_stan/original/2X/b/b6acad3d9f2cd8e7b925c90809f9fb4e4b005a4f.pdf 
 

# Data.frame --------------------------------------------------------
# data.frame with only "rep" column (it's not sort)
indv_gr <- read_tsv("01_RawData/individual_strains_growth_curves_filtered.tsv")
indv_gr <- as.data.frame(indv_gr)

# data.frame with aditional column for replica (ord_replica - replicas sorted)
indv_gro <- read_tsv("01_RawData/modified_individual_strain_growth_curves_ord_replica - individual_strains_growth_curves_filtered.tsv")
indv_gro <- as.data.frame(indv_gro)

# stan_toutputs -- stan model performed with vectors 
outputs_lognormal <- readRDS("03_Output/stan_func_rk_vectors_lognormal")
outputs_normal <- readRDS("03_Output/stan_func_rk_vectors_normal")


# Test distribution ----------------------------------------------------------
CH450_T47 <- indv_gro %>% 
  filter(Cepa == "CH450", # filter each strain 
         temp == 42) %>%  # filter by temp 
  arrange(ord_replica) %>% 
  pull(OD_real)

descdist(CH450_T47, discrete = FALSE)


# Pipeline -------------------------------------------------------------------
# r & k priors from heatmap testing 
stan_rkvectors2 <- stan_ccfunct(df = indv_gro, temp_col = "temp", replica_col = "ord_replica", 
                           strain_col = "Cepa", sample_byh = "hr", interest_col = "OD_real", 
                           inits_list = rk_valz_cc, niterations = 1500, nchains = 2, sigma_val = 0.05)
saveRDS(stan_rkvectors2, file = "03_Output/stan_func_rk_vectors_normal")

# r & k priors (slope and k visualization testing)

rk_stanslope <- stan_rslope(df = indv_gro, xf = rslopedf, temp_col = "temp", replica_col = "ord_replica", 
                            strain_col = "Cepa", sample_byh = "hr", interest_col = "OD_real", niterations = 1500,
                            nchains = 2, sigma_val = 0.05, model = "02_Scripts/loglv_mod.stan", klimit = 2,
                            kcol = "kinit", slope_col = "slope")

saveRDS(rk_stanslope, file ="03_Output/rk_slope")










# ODE check -----------------------------------------------------------------

# COMPARE STAN VS GROWTH CURVES
curvecomparison <- stanvsgw(rk_stanslope, seq(0, 18, 0.1), 0.001, indv_gro, "hr", "Cepa", "temp", "OD_real", 
                            n_samples = 60)
curvecomparison$CH29_T30
curvecomparison$CH111_T37
curvecomparison$CH447_T42


# GROWTH CURVES ONLY 
real_gwtcurves <- growth_curves_func(df = indv_gro, temps = c(30,37,42), strain_col = "Cepa", temp_col = "temp",
                   samplebyh = "hr", interest_column = "OD_real", replica_col = "ord_replica")
real_gwtcurves$CH29_T30
real_gwtcurves$CH29_T37
real_gwtcurves$CH29_T42


