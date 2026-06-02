# Author: Mayra Beatriz Mendoza Velázquez 
# Gompertz model fitting 

# install.packages("MicrobialGrowth")
library(MicrobialGrowth)

# data-------------------------------------------------------------------------
mch2330 <- indv_gro %>%
  filter(Cepa == "CH450", temp == 30) %>%
  arrange(ord_replica, hr) %>% 
  pull(OD_real) %>%
  matrix(ncol = 3, nrow = 10) %>%
  as.data.frame() %>%           
  mutate(time = seq(0,18, by = 2))

# Microbial growth - GOMPERTZ MODEL 
gompertz_ch2330_r1 <- MicrobialGrowth(mch2330$time, mch2330$V1, model = "baranyi")
gompertz_ch2330_r1$f$confint()  # to get confidence interval 

gompertz_ch2330_r2 <- MicrobialGrowth(mch2330$time, mch2330$V2, model = "gompertz")
gompertz_ch2330_r2$f$confint()  # to get confidence interval 

gompertz_ch2330_r3 <- MicrobialGrowth(mch2330$time, mch2330$V3, model = "rosso")
gompertz_ch2330_r3$f$confint()  # to get confidence interval 
?MicrobialGrowth

# plotting 
plot(gompertz_ch2330_r1, coefficients.args = list(cex = 0.75), main = "Regression")
plot(gompertz_ch2330_r2, coefficients.args = list(cex = 0.75), main = "Regression")
plot(gompertz_ch2330_r3, coefficients.args = list(cex = 0.75), main = "Regression")



# Testing all strains #------------------------------------------------------------

gompertz_mgpackage <- function (df, temp_col, replica_col, strain_col, samplebyh, interest_col){
  
  library(tidyverse)
  library(MicrobialGrowth)
  
  spps <- unique(df[[strain_col]])  
  ntemps <- sort(unique(df[[temp_col]]), decreasing = FALSE) 
  
  gomp_out <- list()
  plot_list <- list()
  
  for (m in 1:length(spps)) {
    for (o in 1:length(ntemps)) {
      
      strain_name <- paste0(spps[m], "_T", ntemps[o])
      
      df_complete <- df %>% 
        filter(.data[[strain_col]] == spps[m], 
               .data[[temp_col]] == ntemps[o]) %>%  
        arrange(.data[[replica_col]]) %>%
        select(.data[[samplebyh]], .data[[interest_col]])
      

      fit <- MicrobialGrowth(x = times_vec, y = obs_vec, model = "gompertz")
      
      gomp_out[[strain_name]] <- fit
      
      df_pred <- tibble(
        time    = times_vec,
        OD_real = obs_vec,
        OD_pred = fit$data$predicted
      )
      
      params <- fit$parameters
      

      g <- ggplot(df_pred, aes(x = time)) +
        geom_point(aes(y = OD_real), color = "black", alpha = 0.5, size = 1.8) +
        geom_line(aes(y = OD_pred), color = "dodgerblue3", linewidth = 1.2) +
        theme_classic() +
        labs(
          title = paste("Gompertz Fit:", strain_name),
          subtitle = sprintf("Lag (lambda): %.2f h | Tasa (mu): %.3f h^-1 | Máx (A): %.2f", 
                             params["lambda"], params["mu"], params["A"]),
          x = "Time (hr)", y = "Absorbance (OD600)"
        )
      
      print(g)
    }
  }
  
  return(list(plots = plot_list, gompertz_output = gomp_out))
}


gomp_mgcc <- gompertz_mgpackage(df = indv_gro, temp_col = "temp", replica_col = "ord_replica",
                                strain_col = "Cepa", samplebyh = "hr", interest_col = "OD_real")

plot(gomp_mgcc$plot_list)


x <- indv_gro %>%
  select(hr, OD_real)
class(x)
