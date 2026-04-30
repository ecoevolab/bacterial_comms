# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Functions 

# Stan functions ---------------------------------------------------------------

# 1. Initial values standarized 


# Stan function for logistic LV model fitting 
loglv <- '
functions {
  vector lvfnc(real t,           // time
               vector z,         // state (OD in this particular case)
               real r,           // r
               real k) {         // k
    
    int S = num_elements(z);
    vector[S] dzdt;
    for (j in 1:S) {
      dzdt[j] = r * z[j] * (1 - z[j] / k);
    }
    return dzdt;
  }
}

data {
  int<lower=1> N; 
  int<lower=1> S;
  array[N] real ts;           // time series [without the first time measure]
  vector<lower=0>[S] y0;      // observed initial state 
  array[N, S] real<lower=0> y;      // population measures
}

parameters {
  real<lower=0> r;
  real<lower=0> k;
  vector<lower=0>[S] z0;     // initial state (estimation)
  real<lower=0> sigma;
}

transformed parameters {
  // using vectors instead of real[]
  array[N] vector[S] z = ode_rk45(lvfnc, z0, 0.0, ts, r, k);
}

model {
  // Priors
  r ~ normal(0.3, 0.05);
  k ~ normal(0.8, 0.05);
  sigma ~ lognormal(0, 0.05); # low error 

  // Prior for the initial state based on the values
  z0 ~ normal(y0, 0.1);
}
'

# For standarized initial values 
stan_initial_standarized <- function(S_val){
  list(
    r = 0.3,
    k = 0.8, 
    z0 = as.array(rep(0.001, times = S_val)),
    sigma = 0.1
  )
}

# Stan function with the Log LV model 
stan_ccfunct <- function (df, temp_col, replica_col, strain_col, sample_byh, interest_col, time_series, time_alternative, niterations, nchains){
  
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
  
  
  # Generating initial values and every data.frame for the stan input 
  init_v <- list()
  m_data <- list()
  stan_input <- list()
  ts_vector <- list()
  
  for (q in seq_along(vector_freplica)) {
    
    # extract the first row of every data.frame to get the initial values for stan
    init_v <- as.numeric(vector_freplica[[q]][1, ])
    
    # Extract the rest of the data 
    m_data[[q]] <- vector_freplica[[q]][-1, , drop = FALSE]
    
    # if - to identify specific variations in the data.frame for the time_series 
    
    if (nrow(vector_freplica[[q]]) == length(time_alternative)) {
      
      ts_vector <- time_alternative[-1] 
      
    } else {
      
      ts_vector <- time_series[-1]
      
    }
    
    # That way it doesn't matter if we have another time_series alternative 
    
    # Create the stan input 
    stan_input[[names(vector_freplica)[q]]] <- list(
      N = nrow(m_data[[q]]), 
      S = ncol(m_data[[q]]),
      ts = ts_vector, 
      y0 = init_v,
      y  = m_data[[q]]
    )
  }
  
  stan_output <- list()
  
  for (r in seq_along(stan_input)){
    
    current_S <- stan_input[[r]]$S
    
    init_fun <- function() {
      list(
        r = 0.3,
        k = 0.8, 
        z0 = as.array(rep(0.001, times = current_S)), # it depends on S's value 
        sigma = 0.1
      )
    }
    
    stan_output[[names(stan_input)[r]]] <- stan(
      model_code = loglv, # here it is the stan function i created earlier
      data = stan_input[[r]], # stan_input 
      save_dso = FALSE, 
      iter = niterations,  # iterations
      chains = nchains,   # n. chains 
      init = init_fun) # this function is available in the "Functions" script
    
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
 
