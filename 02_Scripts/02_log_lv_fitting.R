# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Logistic Lotka-Volterra fitting 

library(rstan)
library(readr)
library(dplyr)
library(deSolve)

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
stan_toutputs <- readRDS("bacterial_comms/03_Output/stan_func_rk_vectors")

# stan model performed with matrix 
stan_matrizout <- readRDS("bacterial_comms/03_Output/rk_by_temp_and_strain_withoutCH29")

# Pipeline -------------------------------------------------------------------
pruebafunc <- stan_ccfunct(df = indv_gro, temp_col = "temp", replica_col = "ord_replica", 
                           strain_col = "Cepa", sample_byh = "hr", interest_col = "OD_real", 
                           niterations = 1500, nchains = 3, rvalin = 0.3, kvalin = 0.8, sigma_val = 0.05)

# saveRDS(pruebafunc, file = "03_Output/stan_func_rk_vectors")

# ODE check -----------------------------------------------------------------


# CHECKING ONE BY ONE CODE (STAN VS OD SAMPLING)
# rmean & kmean specific for each list in the outputs list 
# change the following $0CH23_T30 for the replica you want to plot 

rmean <- summary(stan_toutputs$CH23_T42)$summary["r", "mean"]
kmean  <- summary(stan_toutputs$CH23_T42)$summary["k", "mean"]

# parameters 
parameters <- c(
  r = rmean,
  k = kmean
)

# times & initial z value 
timess <- seq(0, 18, 0.1)
init <- c(z = 0.001)


odeCH23 <- ode(y = init, times = timess, func = logistic_eq, parms = parameters)
plot(odeCH23[,1], odeCH23[,2], type = "l", col = "pink", lwd = 2, ylim = c(0, 2),
     main = "Comparación Real: Datos vs Modelo (Bacillus altitudinis - 42°C)")

points(individual_growth$hr[individual_growth$Cepa == "CH23" & individual_growth$temp == 42], 
       individual_growth$OD_real[individual_growth$Cepa == "CH23" & individual_growth$temp == 42], 
       pch = 16, col = "red4")

# USING THE FUNCTION I CREATED TTO COMPARE STAN VS GROWTH CURVES SAMPLES 
curvecomparison <- stanvsgw(stan_toutputs, seq(0, 18, 0.1), 0.001, individual_growth, "hr", "Cepa", "temp", "OD_real")
curvecomparison$CH23_T30

