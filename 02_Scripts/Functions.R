# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Functions 

# Stan functions ---------------------------------------------------------------

# Stan function for logistic LV model fitting 
# Available in the 02_Scripts folder 
# source = https://mc-stan.org/docs/2_26/functions-reference/slicing-and-blocking-functions.html

# For standarized initial values 
# Stan function with the Log LV model 
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

# get real od function ---------------------------------------------------------
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

# Community isolation 

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

# Sparcc for multiple commuinities  - used for testing 
sparcc_inf <- function (list_wcoms, pval){
  
  infernt <- list()
  mt <- "sparcc"
  for (i in seq_along(list_wcoms)){
    infernt[[i]] <- net_inference(taxa_abs = t(list_wcoms[[i]]), method = mt, p = pval)
    
  }
  
  return(infernt)
}
 
