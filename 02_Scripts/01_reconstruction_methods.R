# Author: Mayra Beatriz Mendoza Velazquez
# Title: Exploration of reconstruction methods

library(readr)
library(mlBioNets)
library(dplyr)
library(tibble)
library(SpiecEasi)

# Load tsv archives ----------------------------------------------------------------
metadata <- as.data.frame(read_tsv(file = "bacterial_comms/01_RawData/metadata_clean.tsv"))
metadata[is.na(metadata)] <- 0
f_clean <- as.data.frame(read_tsv("bacterial_comms/01_RawData/f_clean.tsv"))
rzcompositiondata <- read.csv("bacterial_comms/01_RawData/rzcomposition.csv")

## Community isolation -------------------------------------------------------------
comms_rhiz <- unique(metadata$community) # to use all the community values for subsetting the df 
temps_rhiz <- c(28, 32) # without the 0 temperature, because it is already being used in the function 

# Using function 
rz_communities <- community_isolation(rcommunities = comms_rhiz , sample = temps_rhiz , sample_col = "temp", 
                                      rcommunities_col = "community", df = metadata, arrangev = "day", 
                                      interest_column = "label", dfwvals = f_clean, 
                                      composition_df = rzcompositiondata)

    
rz_communities # List with all the communities and their respective composition/abundances 


##### Testing network inference algortihms -----------------------------------------

#### SPARCc ------------------------------------------------------------------------

# 1. Calculate the p value 
sparcc_rhizobial <- sparcc_inf(rz_communities, pval = 0.5)



  