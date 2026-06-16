# Network inference using CVP 
# Author: Mayra Beatriz Mendoza Velazquez 

library(tidyverse)
library(igraph)

# Data -------------------------------------------------------------------------
f_clean <- read_tsv("01_RawData/rhiz_f_clean.tsv")
metadata <- read_tsv("01_RawData/rhiz_metadata_clean.tsv")

f_clean <- f_clean %>%
  rename(Strain = row.names) %>%
  mutate(Strain = replace(Strain, Strain == "NS_042g_27F", "ST00042")) %>%
  mutate(Strain = replace(Strain, Strain == "NS_164C_27F", "ST00164")) %>%
  mutate(Strain = replace(Strain, Strain == "NS_110C_1_27F", "ST00110"))

rzcompositiondata <- read.csv("01_RawData/rzcomposition.csv")
rzcompositiondata <- rzcompositiondata %>%
  mutate(across(everything(), ~ case_when(
    .x == "NS_042g_27F"   ~ "ST00042",
    .x == "NS_164C_27F"   ~ "ST00164",
    .x == "NS_110C_1_27F" ~ "ST00110",
    TRUE ~ .x
  )))

f_clean <- f_clean %>%
  pivot_longer(-Strain,
               names_to = "label",
               values_to = "count")

meta_count <- metadata %>%
  left_join(f_clean %>%
              group_by(label) %>%
              summarise(depth = sum(count),
                        .groups = "drop"),
            by = "label")

f_clean_tocvp <- f_clean %>%
  right_join(meta_count %>% select(community, label, depth), by = "label")





# Nets inference --------------------------------------------------------------- 
n_comms <- unique(f_clean_tocvp$community)

rhiz_complete_tables <- list ()
rhiz_complete_nets <- list()
rzcompositiondata[[1]]


for (i in 1:length(n_comms)){
  
  current_community <- n_comms[i]
  
  current_strains <- rzcompositiondata[[current_community]] %>% 
    na.omit()
  
  
  # get the current community table for lm 
  rhiz_complete_tables[[i]] <- f_clean_tocvp %>%
    filter(community == current_community) %>%
    filter(Strain %in% current_strains) %>%
    pivot_wider(
      names_from = Strain, 
      values_from = count, 
      values_fill = 0
    ) %>%
    select(-depth, -community, -label)
    
  strains_names <- names(rhiz_complete_tables[[i]])
  temporary_net <- list()
  k <- 1

for (Y_cv in  strains_names){
  for(X_cv in  strains_names){
    
    if (Y_cv == X_cv){
      next
    }
    
    # Get the training community for both hypotheses
    all_buty <- setdiff( strains_names, Y_cv)
    all_butxy <- setdiff( strains_names, c(Y_cv, X_cv))
    
    
    # H0: Y = f(Z)
    form_h0 <- as.formula(paste(Y_cv, "~", paste(all_butxy, collapse = " + ")))
    h0_model <- lm(form_h0, data = rhiz_complete_tables[[i]])
    erro_h0 <- sum(residuals(h0_model)^2) # causal strength
    
    
    # H1: Y = f(X,Z)
    form_h1 <- as.formula(paste(Y_cv, "~", paste(all_buty, collapse = " + ")))
    h1_model <- lm(form_h1, data = rhiz_complete_tables[[i]])
    erro_h1 <- sum(residuals(h1_model)^2) # causal_strength
  
    # Causal_strength 
    if(erro_h1 < erro_h0){
      cs <- log(erro_h0/erro_h1)
      newrow <- data.frame(Cause = X_cv, Effect = Y_cv, Causal_Strength = round(cs, 4))
        
      temporary_net[[k]] <- newrow
      k <- k + 1 
      
    } 
    
  }
}

  if (length(temporary_net) > 0) {
    rhiz_complete_nets[[i]] <- bind_rows(temporary_net)
  } else {
    rhiz_complete_nets[[i]] <- data.frame(Cause = character(), Effect = character(), Causal_Strength = numeric())
  }
  
}


# Convert tables into graphs for plotting ----------------------------------------
complete_graphs_cvp <- list()

for (f in 1:length(rhiz_complete_nets)){
 complete_graphs_cvp[[f]] <- graph_from_data_frame(rhiz_complete_nets[[f]], directed = T)
  
}


# plotting -----------------------------------------------------------------------

for (g in 1:length(complete_graphs_cvp)){
  plot(complete_graphs_cvp[[g]])
}
  