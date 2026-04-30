# Author: Mayra Beatriz Mendoza Velazquez 
# Title: Pair-wise network comparison
source("02_Scripts/Functions.R")

library(NetworkComparisonTest)
library(bootnet)


# NCT testing ------------------------------------------------------------

est1 <- estimateNetwork(t(rz_communities[[1]]), default = "IsingFit", tuning = 0)
est2 <- estimateNetwork(t(rz_communities[[2]]), default = "IsingFit", tuning = 0)
set.seed(123)

NCT_b <- NCT(est1, est2, it=10, test.edges=TRUE, 
             edges=list(c(1,2),c(3,5)))

NCT_c <- NCT(est1, est2, paired = TRUE, abs = FALSE, test.edges = TRUE, 
            edges = "all", test.centrality = TRUE, 
            centrality = c("strength"), nodes = "all", it=10)
summary(NCT_c)
summary(NCT_b)
?estimateNetwork
