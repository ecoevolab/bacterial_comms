# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Functions 

# Stan functions ---------------------------------------------------------------

# Stan function with the Log LV model 
# SPECIFIC FOR THE SSE PROCEDURE 
stan_ccfunct <- function (df, temp_col, replica_col, strain_col, sample_byh, interest_col, inits_list, niterations, nchains, sigma_val, model){
  
  # Assigning objects to specific values in the data.frame 
  
  spps <- unique(df[[strain_col]])  # N. of total spps in the data.frame (for the naming of the vectors)
  ntemps <- sort(unique(df[[temp_col]]), decreasing = FALSE) # N. of total temperatures in the data.frame 
  
  # list for the resulting fitting
  stan_output <- list()
  
  # FOR INITIAL R & K AND PRIORS 
  best_inits <- map_df(names(inits_list$grid_output), function(st_name) {
    inits_list$grid_output[[st_name]] %>%
      arrange(sse) %>%   # arrange by best sse
      slice(1) %>%       # get the lowest 
      mutate(id = st_name) %>%
      dplyr::select(id, r_init = r, k_init = k)
  })
  
  for (m in 1:length(spps)) {
    for (o in 1:length(ntemps)) {
      # naming tag 
      strain_name <- paste0(spps[m], "_T", ntemps[o])
      
      # filter the data.frame 
      df_complete <- df %>% 
        filter(.data[[strain_col]] == spps[m], # filter each strain 
               .data[[temp_col]] ==ntemps[o]) %>%  # filter by temp 
        arrange(.data[[replica_col]])     # arrange by the number of replica 
      
      # subset/split the data.frame by replica 
      replica_list <- split(df_complete, df_complete[[replica_col]]) 
      
      # TIME AND OBS VECTORS 
      time_flat <- unlist(lapply(replica_list, function(x) x[[sample_byh]]))
      
      obs_flat <- unlist(lapply(replica_list, function(x) x[[interest_col]]))
      # SPECIFIC PARAMETERS TO THIS FOR LOOP 
      row_id <- best_inits %>%
        filter(id == strain_name)
      
      rvalin <- row_id$r_init
      kvalin <- row_id$k_init
      
        if (kvalin < 1.0) {
          sd_k <- 0.4   
          sd_r <- 0.3
          
        } else {
          sd_k <- 1.0 
          sd_r <- 1.5
          
        }
      
      fixed_inits <- lapply(1:nchains, function(id) {
        list(r = rvalin, k = kvalin)
      })
      
      
      # STAN_DATA 
      stan_data <- list(
        S = length(replica_list), # number of replicas 
        totalobs = length(obs_flat),  # total of observations 
        sizes = as.array(sapply(replica_list, nrow)), # size of each replica 
        y_time = time_flat, 
        y_obs = obs_flat, 
        sigma = sigma_val,
        rlog = rvalin,
        klog = kvalin, 
        sdk = sd_k, 
        sdr = sd_r
      )
      
      
      # Stan function 
      # loglv_mod.stan -- stan function 
      
      message(paste("Stan running for:", strain_name))    
      print(paste("Diagnosis Cepa -> r:", rvalin, "K:", kvalin))
      stan_output[[strain_name]] <- stan(
        file = model,
        data = stan_data, 
        iter = niterations,
        chains = nchains, 
        init = fixed_inits, 
        refresh = 50,
        control = list(adapt_delta = 0.95)
      )
      
    }
  }
  return(stan_output)
}

# STAN WITH SLOPE 

stan_rslope <- function (df, xf, temp_col, replica_col, strain_col, sample_byh, interest_col, niterations, 
                         nchains, sigma_val, model, klimit, kcol, slope_col){
  
  # Assigning objects to specific values in the data.frame 
  
  spps <- unique(df[[strain_col]])  # N. of total spps in the data.frame (for the naming of the vectors)
  ntemps <- sort(unique(df[[temp_col]]), decreasing = FALSE) # N. of total temperatures in the data.frame 
  
  # list for the resulting fitting
  stan_output <- list()
  
  for (m in 1:length(spps)) {
    for (o in 1:length(ntemps)) {
      # naming tag 
      strain_name <- paste0(spps[m], "_T", ntemps[o])
      
      # TIME AND OBS VECTORS 
      # filter the data.frame
      df_complete <- df %>% 
        filter(.data[[strain_col]] == spps[m], # filter each strain 
               .data[[temp_col]] ==ntemps[o]) %>%  # filter by temp 
        arrange(.data[[replica_col]])     # arrange by the number of replica 
      
      # subset/split the data.frame by replica 
      replica_list <- split(df_complete, df_complete[[replica_col]]) 
      
      # TIME AND OBS VECTORS 
      time_flat <- unlist(lapply(replica_list, function(x) x[[sample_byh]]))
      
      obs_flat <- unlist(lapply(replica_list, function(x) x[[interest_col]]))
      
      
      # INITIAL VALUES & K LIMITS FOR STAN TO TEST
      kvalin <- xf %>% 
        filter(.data[[strain_col]] == spps[m], 
               .data[[temp_col]] == ntemps[o]) %>%
        pull(.data[[kcol]]) %>%
        as.numeric()
      
      # sd for stan to test 
      if (kvalin < 1.0) {
        sd_k <- 0.4   
        sd_r <- 1.05
        
      } else {
        sd_k <- 0.4 
        sd_r <- 1.5
        
      }
      
      rslope <- xf %>% 
        filter(.data[[strain_col]] == spps[m], 
               .data[[temp_col]] == ntemps[o]) %>%
        pull(.data[[slope_col]]) %>%
        as.numeric()
      

      fixed_inits <- lapply(1:nchains, function(id) {
        list(
          r = rslope, 
          ka = 0.05,
          yield = kvalin * 1.05)
      })
      
      
      # STAN_DATA 
      stan_data <- list(
        S = length(replica_list), # number of replicas 
        totalobs = length(obs_flat),  # total of observations 
        sizes = as.array(sapply(replica_list, nrow)), # size of each replica 
        
        y_time = time_flat, 
        y_obs = obs_flat, 
        
        sigma = sigma_val,
        rin = rslope,
        kin = kvalin, 
        kfin = klimit, 
        sdk = sd_k, 
        sdr = sd_r
      )
      
      
      # Stan function 
      # loglv_mod.stan -- stan function 
      
      message(paste("Stan running for:", strain_name, "r:", rslope, "K:", kvalin))    

      stan_output[[strain_name]] <- stan(
        file = model,
        data = stan_data, 
        iter = niterations,
        chains = nchains, 
        init = fixed_inits, 
        refresh = 50 
      )
      
    }
  }
  return(stan_output)
}


# logistic ode 

logistic_ode <- function(t, x, parms) {
  r <- parms$r
  k <- parms$k
  
  dxdt <- r * x * (1 - x / k)
  
  return(list(x = dxdt))
}


# compare obs - pred 
compare_obs_pred <- function(r, k , Obs){
  Obs %>%
    map_dfr(function(obs,r,k){
      sim <- ode(y = obs$obs[1],
                 times = obs$hrs,
                 func = logistic_ode,
                 parms = list(r = r, k = k))[,2]
      
      ss <- sum((obs$obs - sim) ^2)
      n <- length(sim)
      
      tibble(ss = ss,
             n = n)
      
    }, r = r, k = k) %>%
    summarise(ss_tot = sum(ss),
              n_tot = sum(n)) %>%
    mutate(sse = ss_tot / n_tot) %>%
    dplyr::select(sse) %>%
    unlist
}


# To get the r slope (lm model)
r_slope <- function(df, xf, strain_col, temp_col, x1_col, x2_col, samplebyh, interest_col) {
  
  # lists 
  y_in_list  <- list()
  y_fin_list <- list()
  slopes_list <- list()
  lm_models   <- list()
  
  # iterations specific for each strain 
  for (i in 1:nrow(xf)) {
    
    current_strain <- xf[[strain_col]][i]
    current_temp   <- xf[[temp_col]][i]
    x1             <- xf[[x1_col]][i]
    x2             <- xf[[x2_col]][i]
    
    
    # get the max value for x1
    y1 <- df %>% 
      filter(.data[[strain_col]] == current_strain,
             .data[[temp_col]] == current_temp,
             .data[[samplebyh]] == x1) %>%
      pull(.data[[interest_col]]) %>% 
      max(na.rm = TRUE)
    
    # get the max value for x2 
    y2 <- df %>% 
      filter(.data[[strain_col]] == current_strain,
             .data[[temp_col]] == current_temp,
             .data[[samplebyh]] == x2) %>%
      pull(.data[[interest_col]]) %>% 
      max(na.rm = TRUE)
    
    # save the y1 & y2 vals 
    y_in_list[[i]]  <- y1
    y_fin_list[[i]] <- y2
    
    # lineal model for eachs strain
    if (!is.na(y1) && !is.na(y2)) {
      x_vals <- c(x1, x2)
      y_vals <- c(y1, y2)
      
      fit <- lm(y_vals ~ x_vals)
      
      lm_models[[i]]   <- fit
      slopes_list[[i]]  <- coef(fit)[2]
      
    }}
  
  # new columns to include in the data.frame 
  xf[["yin"]]    <- y_in_list
  xf[["yfin"]]   <- y_fin_list
  xf[["slope"]]   <- slopes_list
  
  return(xf)
}

# Obs filtering for r & k prior estimation 
rk_prior_testing <- function (df, temp_col, replica_col, strain_col, sample_byh, interest_col){
  
  # Assigning objects to specific values in the data.frame 
  
  spps <- unique(df[[strain_col]])  # total n. of spps in the data.frame (for the naming of the vectors)
  ntemps <- sort(unique(df[[temp_col]]), decreasing = FALSE) # N. of total temperatures in the data.frame 
  
  # list for the resulting fitting
  obs_filtered <- list()
  grid_results <- list()
  
  for (m in 1:length(spps)) {
    for (o in 1:length(ntemps)) {

        # filter the data.frame 
        df_complete <- df %>% 
          filter(.data[[strain_col]] == spps[m], # filter each strain 
                 .data[[temp_col]] ==ntemps[o]) %>%  # filter by temp 
          arrange(.data[[replica_col]], .data[[sample_byh]])     # arrange by the number of replica 
        
        strain_temptag <- paste0(spps[m], "_T", ntemps[o])
        replica_list <- df_complete %>%
                        group_split(.data[[replica_col]]) %>%
                        map(~{
                           .x %>%
                           dplyr::select(hrs = all_of(sample_byh), obs = all_of(interest_col))
                           })
        obs_filtered[[strain_temptag]] <- replica_list
        
        Res <- expand_grid(r = seq(0.1, 2, by = 0.1),
                           k = seq(0.5, 2, by = 0.1))
        
        # obs == replica_list 
        Res$sse <- Res %>% 
          dplyr::select(r, k) %>%
          pmap_dbl(.f = compare_obs_pred, Obs = replica_list)
        
        # keep the results 
        grid_results[[strain_temptag]] <- Res %>% arrange(sse)
        
        # plotting 
        p <- ggplot(Res, aes(x = r, y = k)) +
          geom_tile(aes(fill = sse)) +
          scale_fill_gradient2(low = "blue", mid = "lightyellow", high = "red",
                               midpoint = mean(Res$sse)) +
          labs(title = strain_temptag) +
          theme_classic()
        
        print(p) 
        
   }
  }
  return(list(obs = obs_filtered, grid_output = grid_results))
}

# logistic growth (deSolve) / for testing stan model
logistic_eq <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dz <- r * z * (1 - z / k)
    return(list(c(dz)))
  })
}


# Growth curves ---------------------------------------------------------------

# Get each growth curve (cuatro cienegas)
growth_curves_func <- function(df, temps, strain_col, temp_col, samplebyh, interest_column, replica_col){
  
  strains <- unique(df[[strain_col]])
  all_curvesbyt <- list()
  tempt <- as.numeric(temps)
  
  for (i in 1:length(strains)){
    for (f in 1:length(tempt)){
      title <- paste0(strains[i], "_T", tempt[f])
      
      all_curvesbyt[[title]] <- df[df[[strain_col]] == strains[i] & df[[temp_col]] == tempt[f], ]
    }
  }
  
  plot_list <- list()
  
  for (g in 1:length(all_curvesbyt)){
    current_name <- names(all_curvesbyt)[g]
    gcurveplot <- ggplot(all_curvesbyt[[g]], aes(x = .data[[samplebyh]], 
                                                 y = .data[[interest_column]], 
                                                 group = .data[[replica_col]], 
                                                 colour = as.factor(.data[[replica_col]]))) +
      geom_point() + 
      geom_line(na.rm = TRUE) + 
      theme_classic() +
      xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
      labs(title = paste("Growth curve:", current_name))
    
    plot_list[[current_name]] <- gcurveplot
  }
  return(plot_list)
}


# Compare stan results vs the actual growth curves 
stanvsgw <- function(outs, timesf, initialv, realdata, samplebyh, strain_col, temp_col, interest_col, n_samples){
  
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(tidyr)
  
  plot_list <- list()
  names_outs <- names(outs) 
  
  for (f in 1:length(outs)){
    
    current_tag <- names_outs[f]
    rkvals <- as.data.frame(outs[[f]], pars = c("r", "k"))
    
    # select an x number of samples to graph/plot 
    set.seed(123) 
    vz <- sample(1:nrow(rkvals), min(n_samples, nrow(rkvals)))
    samples_df <- rkvals[vz, ]
    
    # create curves 
    n0 <- initialv 
    curves <- map_df(1:nrow(samples_df), function(i) {
      r_i <- samples_df$r[i]
      k_i <- samples_df$k[i]
      
      # get the logistic curve 
      tibble(
        time = timesf,
        OD_sim = k_i / (1 + ((k_i - n0) / n0) * exp(-r_i * timesf)),
        iter = i
      )
    })
    
    # mean line (just to see it and compare with the rest 
    r_mean <- mean(rkvals$r)
    k_mean <- mean(rkvals$k)
    
    meanline <- tibble(
      time = timesf,
      OD_sim = k_mean / (1 + ((k_mean - n0) / n0) * exp(-r_mean * timesf))
    )
    
    # filter real data 
    split_name <- strsplit(current_tag, split = "_")[[1]]
    strainm <- split_name[1]
    temps   <- as.numeric(gsub("T", "", split_name[2]))
    
    real_subset <- realdata[realdata[[strain_col]] == strainm & realdata[[temp_col]] == temps, ]
    
    # 5. Plotting
    g <- ggplot() +
      
      # group of lines 
      geom_line(data = curves, aes(x = time, y = OD_sim, group = iter), 
                color = "deeppink4", alpha = 0.1) +
      
      # mean line - higher linewidth so we can see it better 
      geom_line(data = meanline, aes(x = time, y = OD_sim), 
                color = "deeppink", linewidth = 1.2) +
      
      # real data points 
      geom_point(data = real_subset, aes(x = .data[[samplebyh]], y = .data[[interest_col]]), 
                 color = "red4", size = 2, alpha = 0.7) +
      
      theme_classic() +
      labs(title = paste("Posterior Predictive Check:", current_tag),
           subtitle = sprintf("Media: r = %.3f | k = %.3f", r_mean, k_mean),
           x = "Time (hr)", y = "Absorbance (OD600)")
    
    plot_list[[current_tag]] <- g
  }
  
  return(plot_list)
}


# Creates a list of matrix with every OD_real value (separated by strain and temperature)
get_realod <- function (df, temp_col, replica_col, strain_col, sample_byh, interest_col){
  
  # assigning objects to specific values in the data.frame 
  
  spps <- unique(df[[strain_col]]) 
  ntemps <- sort(unique(df[[temp_col]]), decreasing = FALSE)
  
  ntemps_numeric <- as.numeric(ntemps) 
  ntemps_character <- as.character(ntemps)
  
  vector_freplica <- list() 
  p <- 1
  
  
  for (m in 1:length(spps)) {
    for (o in 1:length(ntemps_numeric)) {
      
      df_complete <- df[df[[strain_col]] == spps[m] & df[[temp_col]] == ntemps_numeric[o], ]
      
      timepoints <- length(unique(df_complete[[sample_byh]]))
      ntotalc <- nrow(df_complete) / timepoints
      
      if (nrow(df_complete) > 0){
        
        df_filtered <- df_complete %>% 
          arrange(.data[[replica_col]], .data[[sample_byh]]) %>%
          pull(.data[[interest_col]])
        
        df_matrix <- matrix(df_filtered, ncol = ntotalc)
        colnames(df_matrix) <- rep(ntemps_numeric[o], times = ntotalc)
        
        vector_freplica[[p]] = df_matrix
        
        names(vector_freplica)[p] <- paste0(spps[m], "_T", ntemps_numeric[o])
        p <- p + 1
      }
    }
  }
  
  return(vector_freplica)
  
}
# Reconstruction methods -------------------------------------------------------

# Get each community abundances (separated by temperature) 

community_isolation <- function(rcommunities, sample, sample_col, rcommunities_col, df, arrangev, interest_column, dfwvals, composition_df){
  
  # create an empty list   
  community_list <- list()
  k <- 1
  
  # subset the commuinity values based on the day and community name 
  
  for (i in 1:length(rcommunities)){
    for (x in 1:length(sample)){
      
      community_list[[k]] <- subset(df, df[[sample_col]] %in% c(0, sample[x])  & df[[rcommunities_col]] == rcommunities[i]) %>% 
        arrange(.data[[arrangev]]) %>% 
        pull(.data[[interest_column]]) %>% 
        as.character()
      
      k <- k + 1
    }
  }
  
  # create an empty list for the abundances 
  abundances_tables <- list()
  
  # Based on the column names, select the values from the data.frame with sample measures  
  
  for (id in seq_along(community_list)) {
    abundances_tables[[id]] <- dfwvals %>%
      dplyr::select(1 | all_of(community_list[[id]])) %>%
      column_to_rownames(var = "row.names") # Asumming the first column belongs to the row names 
    
  }
  
  abnds_filtered <- list()
  
  for (g in seq_len(ncol(composition_df))) {
    
    # select the community name / for arranging the list 
    comm_name <- colnames(composition_df)[g]
    
    # select the specific bacterial id for the selected community 
    spp <- composition_df[, g]
    spln <- 1:length(spp)
    
    
    # indexes [to select in the list of data.frames]
    idx <- (2*g - 1):(2*g)
    
    for (j in seq_along(idx)) {
      
      k <- idx[j]
      
      abnds_filtered[[paste0(comm_name, "_", sample[j])]] <- abundances_tables[[k]] %>% 
        filter(row.names(abundances_tables[[k]]) %in% spp)
      
      
    }
  }
  
  return(abnds_filtered)
  
}

# Sparcc for multiple commuinities  - used for algorithm testing 
sparcc_inf <- function (list_wcoms, pval){
  
  infernt <- list()
  mt <- "sparcc"
  for (i in seq_along(list_wcoms)){
    infernt[[i]] <- net_inference(taxa_abs = t(list_wcoms[[i]]), method = mt, p = pval)
    
  }
  
  return(infernt)
}
 
