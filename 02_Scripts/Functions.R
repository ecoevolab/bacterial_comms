# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Functions 

# Stan functions ---------------------------------------------------------------

# 1. Stan function with the Log LV model 
stan_ccfunct <- function (df, temp_col, replica_col, strain_col, sample_byh, interest_col, niterations, nchains, rvalin, kvalin, sigma_val){
  
  # Assigning objects to specific values in the data.frame 
  
  spps <- unique(df[[strain_col]])  # N. of total spps in the data.frame (for the naming of the vectors)
  ntemps <- sort(unique(df[[temp_col]]), decreasing = FALSE) # N. of total temperatures in the data.frame 
  
  # list for the resulting fitting
  stan_output <- list()
  
  for (m in 1:length(spps)) {
    for (o in 1:length(ntemps)) {
      
      # filter the data.frame 
      df_complete <- df %>% 
        filter(.data[[strain_col]] == spps[m], # filter each strain 
               .data[[temp_col]] ==ntemps[o]) %>%  # filter by temp 
        arrange(.data[[replica_col]])     # arrange by the number of replica 
      
       # subset/split the data.frame by replica 
      replica_list <- split(df_complete, df_complete[[replica_col]]) 
      
      # create the time and obs vectors 
      time_flat <- unlist(lapply(replica_list, function(x) x[[sample_byh]]))
      
      obs_flat <- unlist(lapply(replica_list, function(x) x[[interest_col]])) 


      # generate the stan_data that's going to be used in the stan function 
      stan_data <- list(
        S = length(replica_list), # number of replicas 
        totalobs = length(obs_flat),  # total of observations 
        sizes = as.array(sapply(replica_list, nrow)), # size of each replica 
        y_time = time_flat, 
        y_obs = obs_flat, 
        sigma = sigma_val
      )
      
      # assign a name to the object 
      strain_name <- paste0(spps[m], "_T", ntemps[o])
      
      #Initial values for the stan model 
      init_fun <- function(){
        list(
        r = as.numeric(rvalin), 
        k = as.numeric(kvalin))
      }
      
      # Stan function 
      # loglv_mod.stan -- stan function 
      
      stan_output[[strain_name]] <- stan(
        file = "02_Scripts/loglv_mod.stan",
        data = stan_data, 
        iter = niterations,
        chains = nchains, 
        init = init_fun, 
        refresh = 50 
      )
      
    }
  }
  return(stan_output)
}


# 2. logistic growth (deSolve) / for testing stan model
logistic_eq <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dz <- r * z * (1 - z / k)
    return(list(c(dz)))
  })
}


# Growth curves ---------------------------------------------------------------

# 3. Get each growth curve (cuatro cienegas)
growth_curves_func <- function(df, temps, strain_col, temp_col, samplebyh, interest_column, replica_col){
  
  strains <- sort(unique(df[[strain_col]]), decreasing = F)
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


# 4. Compare stan results vs the actual growth curves 
stanvsgw <- function(outs, timesf, initialv, realdata, samplebyh, strain_col, temp_col, interest_col){
  
  library(ggplot2)
  library(deSolve)
  
  init <- c(z = initialv)
  plot_list <- list()
  names_outs <- names(outs) 
  
  for (f in 1:length(outs)){
    
    # Extract "outs" summary (to get the specific r & k values)
    current_tag <- names_outs[f]  # get the 1..2..3.. outs names
    summ <- summary(outs[[f]])$summary # get the summary for the specific name 
    
    # get r & k values 
    r_val <- as.numeric(summ["r", "mean"])
    k_val <- as.numeric(summ["k", "mean"])
    p_actual <- c(r = r_val, k = k_val) # to use in the deSolve equation
    
    # solve ode & change column names in the ode object 
    outsolve <- as.data.frame(ode(y = init, times = timesf, func = logistic_eq, parms = p_actual))
    colnames(outsolve) <- c("time", "OD_sim")
    
    # Filter real data (actual OD sampling)
    split_name <- strsplit(current_tag, split = "_")[[1]]
    strain_n <- split_name[1]
    temp_v   <- as.numeric(gsub("T", "", split_name[2]))
    
    # this way i get the specific strain & temp values to compare to the ode results 
    real_subset <- realdata[realdata[[strain_col]] == strain_n & realdata[[temp_col]] == temp_v, ]
    
    # Plotting 
    # We can use r_val & k_val because these are numeric values 
    g <- ggplot() +
      geom_line(data = outsolve, aes(x = time, y = OD_sim), color = "deeppink", linewidth = 1) +
      
      geom_point(data = real_subset, aes(x = .data[[samplebyh]], y = .data[[interest_col]]), 
                 color = "red4", alpha = 0.6) +
      
      theme_classic() +
      ylim(0, 1.8) +
      labs(title = paste("Stan vs Real:", current_tag),
           subtitle = sprintf("r = %.3f | k = %.3f", r_val, k_val),
           x = "Time (hr)", y = "Absorbance (OD600)")
    
    plot_list[[current_tag]] <- g
  }
  
  return(plot_list)
}


# 5. Creates a list of matrix with every OD_real value (separated by strain and temperature)
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

# 6. Get each community abundances (separated by temperature) 

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

# 7. Sparcc for multiple commuinities  - used for algorithm testing 
sparcc_inf <- function (list_wcoms, pval){
  
  infernt <- list()
  mt <- "sparcc"
  for (i in seq_along(list_wcoms)){
    infernt[[i]] <- net_inference(taxa_abs = t(list_wcoms[[i]]), method = mt, p = pval)
    
  }
  
  return(infernt)
}
 
