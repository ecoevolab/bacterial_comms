# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Stan model fitting 

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







# Pipeline -------------------------------------------------------------------
# r & k priors from heatmap testing (stan_ccfunct)
stan_rkvectors2 <- stan_ccfunct(df = indv_gro, temp_col = "temp", replica_col = "ord_replica", 
                           strain_col = "Cepa", sample_byh = "hr", interest_col = "OD_real", 
                           inits_list = rk_valz_cc, niterations = 1500, nchains =4 , sigma_val = 0.05)



# r & k priors (slope and k visualization testing) (stan_rslope)

rk_slope_priors <- readRDS("03_Output/rkpriors_slope")
rkmonod_stanslope <- stan_rslope(df = indv_gro, xf = rk_slope_priors, temp_col = "temp", replica_col = "ord_replica", 
                            strain_col = "Cepa", sample_byh = "hr", interest_col = "OD_real", niterations = 1500,
                            nchains = 2, sigma_val = 0.05, model = "02_Scripts/monod_mod.stan", klimit = 2,
                            kcol = "kinit", slope_col = "slope")




# MONOD #
prior_heatmap <- readRDS("03_Output/rkpriors_heatmap")
monod_mdl <-  stan_monodf(df = indv_gro, temp_col = "temp", replica_col = "ord_replica", 
               strain_col = "Cepa", sample_byh = "hr", interest_col = "OD_real", 
               inits_list = prior_heatmap, niterations = 1500, nchains = 3 , sigma_val = 0.05, 
               model = "02_Scripts/monod_mod.stan")
saveRDS(monod_mdl, "03_Output/stan_rk_vectors_heatmap_monod")






# ODE check -----------------------------------------------------------------
monod_output <- readRDS("03_Output/stan_rk_vectors_slope_monod")

# COMPARE STAN VS GROWTH CURVES -- MONOD 
monod_results1 <- monodvsgw(outs = monod_output, timesf = seq(0, 18, by = 0.1), realdata = indv_gro, samplebyh = "hr",
                           strain_col   = "Cepa", temp_col     = "temp", interest_col = "OD_real", n_samples    = 30)
            
monod_results1$CH149a_T42

monod_results2 <- monodvsgw(outs = monod_mdl, timesf = seq(0, 18, by = 0.1), realdata = indv_gro, samplebyh = "hr",
                           strain_col   = "Cepa", temp_col     = "temp", interest_col = "OD_real", n_samples    = 30)

monod_results2






# GROWTH CURVES ONLY 
real_gwtcurves <- growth_curves_func(df = indv_gro, temps = c(30,37,42), strain_col = "Cepa", temp_col = "temp",
                   samplebyh = "hr", interest_column = "OD_real", replica_col = "ord_replica")




















monodvsgw <- function(outs, timesf, realdata, samplebyh, strain_col, temp_col, interest_col, n_samples){
  
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(deSolve) 
  
  
  monod_equations <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      X <- max(X, 1e-6)
      S <- max(S, 1e-6)
      
      mu <- mumax * (S / (ks + S))
      
      dXdt <- mu * X
      dSdt <- - (mu * X) / yield
      
      return(list(c(dXdt, dSdt)))
    })
  }
  
  plot_list <- list()
  names_outs <- names(outs) 
  
  for (f in 1:length(outs)){
    
    current_tag <- names_outs[f]
    
    split_name <- strsplit(current_tag, split = "_")[[1]]
    strainm <- split_name[1]
    temps   <- as.numeric(gsub("T", "", split_name[2]))
    
    real_subset <- realdata %>% 
      filter(.data[[strain_col]] == strainm, .data[[temp_col]] == temps)
    
    if(nrow(real_subset) == 0) next
    
    X0_real <- min(real_subset[[interest_col]])
    if(X0_real <= 0) X0_real <- 0.001
    
    posterior_samples <- as.data.frame(rstan::extract(outs[[f]], pars = c("mumax", "ks", "yield")))
    
    set.seed(123) 
    vz <- sample(1:nrow(posterior_samples), min(n_samples, nrow(posterior_samples)))
    samples_df <- posterior_samples[vz, ]
    
    curves <- map_df(1:nrow(samples_df), function(i) {
      
      state_init <- c(X = X0_real, S = 1.0)
      
      pars <- c(
        mumax = samples_df$mumax[i],
        ks    = samples_df$ks[i],
        yield = samples_df$yield[i]
      )
      
      ode_output <- as.data.frame(ode(y = state_init, times = timesf, func = monod_equations, parms = pars))
      
      tibble(
        time = ode_output$time,
        OD_sim = ode_output$X, 
        iter = i
      )
    })
    
    mumax_mean <- mean(posterior_samples$mumax)
    ks_mean    <- mean(posterior_samples$ks)
    yield_mean <- mean(posterior_samples$yield)
    
    mean_ode <- as.data.frame(ode(
      y = c(X = X0_real, S = 1.0), 
      times = timesf, 
      func = monod_equations, 
      parms = c(mumax = mumax_mean, ks = ks_mean, yield = yield_mean)
    ))
    
    meanline <- tibble(
      time = mean_ode$time,
      OD_sim = mean_ode$X
    )
    
    g <- ggplot() +
      geom_line(data = curves, aes(x = time, y = OD_sim, group = iter), 
                color = "deeppink4", alpha = 0.15) +
      
      geom_line(data = meanline, aes(x = time, y = OD_sim), 
                color = "deeppink", linewidth = 1.2) +
      
      geom_point(data = real_subset, aes(x = .data[[samplebyh]], y = .data[[interest_col]]), 
                 color = "black", size = 1.8, alpha = 0.6) +
      
      theme_classic() +
      labs(
        title = paste("Monod Predictive Check:", current_tag),
        subtitle = sprintf("Medias: mu_max = %.2f h^-1 | Ks = %.3f | Yield = %.2f", 
                           mumax_mean, ks_mean, yield_mean),
        x = "Time (hr)", 
        y = "Absorbance (OD600)"
      )
    
    plot_list[[current_tag]] <- g
  }
  
  return(plot_list)
}
