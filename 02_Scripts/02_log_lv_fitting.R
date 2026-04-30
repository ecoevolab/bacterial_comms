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

# Table interpretation notes --------------------------------------------------
# r & k are the ones that i'm interested in 
# r - fitness 
# k - carrying capacity 
# mean 
# se_mean 
# sd - standard deviation 
# % -> quantiles, confidence intervals ?? these provide an estimate of range where 
# the true value is likely to be 
# n_eff -> number of effective samples / if it is too low the aproximations 
# may lack precision, we want this number to be specially high 
# Rhat - Do the mcmc chains get to the same value? - 1.01 is a good indicator
# if this values goes above from 1.01 or 1.05 then that means the value of interest (r, k or the initial one)
# are not realiable values and we cannot use them in the model 
# Values close to 1.0 indicate that the chains have successfully converged to the same distribution

# Pipeline -------------------------------------------------------------------
# we're going to use the data.frame without CH29 because the sampling is different from the other spps 

indv_withoutCH29 <- indv_gro %>% filter(Cepa != "CH29")
pruebafunc <- stan_ccfunct(df = indv_withoutCH29, temp_col = "temp", replica_col = "rep", 
                           strain_col = "Cepa", sample_byh = "hr", interest_col = "OD_real", 
                           time_series = seq(0, 18, by = 2), time_alternative = c(0, 10, 12, 14, 16, 18), 
                           niterations = 3000, nchains = 4)

rk_valslist_wtCH29 <- saveRDS(pruebafunc, file = "03_Output/rstan_logistic_fitting/rk_withoutCH29_yrepincluded")


# CH29 - 10 times specific -------------------------------------------------------------
time <- seq(0, 18, by = 2)

# OD
rep1_CH29 <- subset(indv_gro, Cepa == "CH29" &
                      ((rep == 1 & temp == 42))) %>%
  arrange(temp) %>%
  pull(OD_real)

rep2_CH29 <- subset(indv_gro, Cepa == "CH23" &
                      ((rep == 3 & temp == 42))) %>%
  arrange(temp) %>%
  pull(OD_real)

# OD [ column are the different replicas]
OD_CH29 <- cbind(
  # replica 1
  rep1_CH29,
  # replica 2
  rep2_CH29
)

# initial values 
in_val <- c(rep1_CH29[1], rep2_CH29[1])
N_val <- length(time)

log_CH29_df <- list(
  
  N = 9,
  S = ncol(OD_CH29),
  ts = time[-1],
  y0 = in_val,
  y = OD_CH29[-1,]
)

S_29_10t <- log_CH29_df$S

init_fun29_10t <- function() {
  list(
    r = 0.3,
    k = 0.8, 
    z0 = as.array(rep(0.001, times = 2)), 
    sigma = 0.1
  )
}

log_CH29_10tfit <- stan(model_code = loglv, 
                    data = log_CH29_df, 
                    save_dso = FALSE, 
                    iter = 3000,  # iterations
                    chains = 4,   # n. chains 
                    init = init_fun29_10t ) # the initial values i want stan to start from 
print(log_CH29_10tfit)
saveRDS(log_CH29_10tfit, file = "03_Output/CH29_T42_10times")

# CH29 - 6 times specific -------------------------------------------------------------
time6t_CH29 <- c(0,10,12,14,16,18)

# OD
rep1_CH29_6t <- subset(indv_gro, Cepa == "CH29" &
                      ((rep == 4 & temp == 42))) %>%
  arrange(temp) %>%
  pull(OD_real)

rep2_CH29_6t <- subset(indv_gro, Cepa == "CH29" &
                      ((rep == 5 & temp == 42))) %>%
  arrange(temp) %>%
  pull(OD_real)

# OD [ column are the different replicas]
OD_CH29_6t <- cbind(
  # replica 1
  rep1_CH29_6t,
  # replica 2
  rep2_CH29_6t
)

# initial values 
in_val_6t <- c(rep1_CH29_6t[1], rep2_CH29_6t[1])
N_val_6t <- length(time6t_CH29)

log_CH29_df_6t <- list(
  
  N = 5,
  S = ncol(OD_CH29_6t),
  ts = time6t_CH29[-1],
  y0 = in_val_6t,
  y = OD_CH29_6t[-1,]
)

S_29_6t <- log_CH29_df_6t$S

init_fun29_6t <- function() {
  list(
    r = 0.3,
    k = 0.8, 
    z0 = as.array(rep(0.001, times = S_29_6t)), 
    sigma = 0.1
  )
}

log_CH29_6fit <- stan(model_code = loglv, 
                    data = log_CH29_df_6t, 
                    save_dso = FALSE, 
                    iter = 3000,  # iterations
                    chains = 4,   # n. chains 
                    init = init_fun29_6t ) # the initial values i want stan to start from 
print(log_CH29_6fit)
saveRDS(log_CH29_6fit, file = "03_Output/CH29_T42_6times")
