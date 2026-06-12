# Cross Validation by community only
# Author: Mayra Mendoza 

library(tidyverse)

# Data.frame abundances 
f_clean <- read_tsv("01_RawData/rhiz_f_clean.tsv")
f_clean <- f_clean %>%
  rename(Strain = row.names) %>%
  mutate(Strain = replace(Strain, Strain == "NS_042g_27F", "ST00042")) %>%
  mutate(Strain = replace(Strain, Strain == "NS_164C_27F", "ST00164")) %>%
  mutate(Strain = replace(Strain, Strain == "NS_110C_1_27F", "ST00110"))
f_clean

# data.frame replica label
metadata <- read_tsv("01_RawData/rhiz_metadata_clean.tsv")
metadata <- metadata %>%
  mutate(sample = str_remove(label, "[AB]$"))
metadata


# data.frame OD values 
syncom_growth <-read_tsv("pilot_syncom_growth_curves.tsv")
syncom_growth <- syncom_growth %>%
  select(Community, temp, OD600, hrs = total_time_h, exp = batch)
syncom_growth


# data.frame with specific community structure 
rzcompositiondata <- read_csv("01_RawData/rzcomposition.csv") 


freqs <- f_clean %>%
  pivot_longer(-Strain,
               names_to = "label",
               values_to = "count")

# generate the "count" of the strain in the media 
meta_count <- metadata %>%
  left_join(freqs %>%
              group_by(label) %>%
              summarise(depth = sum(count),
                        .groups = "drop"),
            by = "label")


freqs <- freqs %>%
  right_join(meta_count %>% select(community, label, depth), by = "label") %>%
  mutate(freq = count / depth) %>%
  select(-depth)

meta_count <- meta_count %>%
  select(-day, -techrep) %>%
  left_join(syncom_growth, by = c("Community", "hrs", "exp", "temp"))


freqsss <- freqs %>%
  right_join(indiv_growth %>% select(OD600), by = "Strain")

Freqs <- Freqs %>%
  left_join(Meta %>% select(label, OD600)) %>%
  mutate(abs1 = freq * OD600) %>%
  select(-OD600)


# Tibble for the CV

CVP_table <- freqs %>%
  select(Strain, label, community, freq) %>% 
  pivot_wider(
    names_from = community, 
    values_from = freq, 
    values_fill = 0
  )


# Cross validation 
commnames <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12")


community_CVP <- list()

for (Y_cv in commnames){
  for(X_cv in commnames){
    
    if (Y_cv == X_cv){
      next
    }
    
    # Get the training community for both hypotheses
    all_buty <- setdiff(commnames, Y_cv)
    all_butxy <- setdiff(commnames, c(Y_cv, X_cv))
    
    
    # H0: Y = f(Z)
    form_h0 <- as.formula(paste(Y_cv, "~", paste(all_butxy, collapse = " + ")))
    h0_model <- lm(form_h0, data = CVP_table)
    erro_h0 <- sum(residuals(h0_model)^2) # causal strength
    
    
    # H1: Y = f(X,Z)
    form_h1 <- as.formula(paste(Y_cv, "~", paste(all_buty, collapse = " + ")))
    h1_model <- lm(form_h1, data = CVP_table)
    erro_h1 <- sum(residuals(h1_model)^2) # causal_strength
    
    erro_j <- erro_h1 - erro_h0
    if(erro_j < 0){
      
    }
    
    # Causal_strength 
    if(erro_h1 < erro_h0){
      cs <- log(erro_h0/erro_h1)
      newrow <- data.frame(Cause = X_cv, Effect = Y_cv, Causal_Strength = round(cs, 4))
        
      community_CVP <- append(community_CVP, list(newrow))
      
    } 
    
  }
}

CVP_ftab <- bind_rows(community_CVP)
CVP_ftab
