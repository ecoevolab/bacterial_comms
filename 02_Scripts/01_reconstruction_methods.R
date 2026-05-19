# Author: Mayra Beatriz Mendoza Velazquez
# Title: Exploration of reconstruction methods

library(readr)
library(mlBioNets)
library(dplyr)
library(tibble)
library(SpiecEasi)
library(textshape)
install.packages("WGCNA")

# Load tsv archives ----------------------------------------------------------------
metadata <- as.data.frame(read_tsv(file = "01_RawData/metadata_clean.tsv"))
metadata[is.na(metadata)] <- 0
f_clean <- as.data.frame(read_tsv("01_RawData/f_clean.tsv"))
rzcompositiondata <- read.csv("01_RawData/rzcomposition.csv")

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

# 1. Isolate only by community (not temperature)
community_list <- list()
k <- 1
# subset the commuinity values based on the day and community name 
for (i in 1:length(comms_rhiz)){

    community_list[[k]] <- subset(metadata, community %in% comms_rhiz[i])  %>% 
      arrange("day") %>% 
      pull("label") %>% 
      as.character()
    
    k <- k + 1
  
}

# create an empty list for the abundances 
abundances_tables <- list()

# Based on the column names, select the values from the data.frame with sample measures  
for (id in seq_along(community_list)) {
  abundances_tables[[id]] <- f_clean %>%
    dplyr::select(1 | all_of(community_list[[id]]))  %>%
    column_to_rownames("row.names") # Asumming the first column belongs to the row names 
  
}

abnds_filtered <- list()

for (g in seq_len(ncol(rzcompositiondata))) {
  
  # select the community name / for arranging the list 
  comm_name <- colnames(rzcompositiondata)[g]
  
  # select the specific bacterial id for the selected community 
  spp <- rzcompositiondata[, g]
  spln <- 1:length(spp)
    
    abnds_filtered[[paste0(comm_name)]] <- abundances_tables[[g]] %>% 
      filter(row.names(abundances_tables[[g]]) %in% spp)
    
}

#### ALGORITHMS #### 

# MRNET # 
# Based on mutual information 

# WGCNA #


#### SPARCc ------------------------------------------------------------------------
library (mlBioNets)
# Undirected
# Based on correlations 
sparcc_rhizobial <- sparcc_inf(abnds_filtered, pval = 0.05)
# the function gets the p and cor matrixes by itself, pval corresponds only to the significance 

#### GENIE3 -----------------------------------------------------------------------
# Directed & weighted 
# Based on 

library(GENIE3)

# weighted matrix 
weightMat <- GENIE3(exprMatr)

# get regulatory links 
linkList <- getLinkList(weightMat)
top_edges <- linkList[1:500, ]
grn_graph <- graph_from_data_frame(d = top_edges, directed = TRUE)

# assigning vertex & edge attributes 
V(grn_graph)$size <- 5
V(grn_graph)$label <- V(grn_graph)$name
V(grn_graph)$label.cex <- 0.6
V(grn_graph)$color <- "lightblue"

# Plot the directed network
plot(grn_graph, 
     layout = layout_with_fr(grn_graph), # Fruchterman-Reingold layout
     vertex.label.dist = 0.5,
     edge.arrow.size = 0.3,
     main = "Inferred Gene Regulatory Network")


  