# Author: Mayra Beatriz Mendoza Velazquez 
# Title: Rstan 

library (bayesplot)
library(ggplot2)
library(hexbin)
library(dplyr)
library(tibble)
library(posterior)
library(tidyr)
library(stringr)
library(rstan)
library(gridExtra)
library(rstanemax)

# Load the data (stan outputs) ---------------------------------------------------
stan_toutputs <- readRDS("bacterial_comms/03_Output/stan_func_rk_vectors")

# Merge two lists (for posterior manipulation) -------------------------------------

# TRACE GRAPHS -------------------------------------------------------------------
# Ref: https://mc-stan.org/rstan/reference/stanfit-method-plot.html

# mcmc_trace ----------------------------------------------------------------
pdf("03_Output/Traceplots_CCStrains.pdf", width = 15, height = 7)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]
  
  pst_cp <- as.array(stan_toutputs[[i]])
  
  color_scheme_set("mix-brightblue-gray")
  
  # adding title to each plot 
  f <- mcmc_trace(pst_cp, pars = c("r", "k", "sigma")) +
    labs(title = paste("Traceplot:", strain),
         subtitle = "Sampling (1500)")
  
  print(f)
}

dev.off()


# plot()----------------------------------------------
pdf("bacterial_comms_output/traceplot_rk_5_5_26.pdf", width = 15, height = 7)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]

  # adding title to each plot 
  x <- plot(stan_toutputs[[i]], plotfun = "trace", pars = c("r", "k"), inc_warmup = TRUE) +
    labs(title = paste("Traceplot:", strain),
         subtitle = "Shaded area: Warmup / Unshaded area: Sampling")
  
  print(x)
}

dev.off()


# PAIRS GRAPHS -------------------------------------------------------------------

# pairs func -------------------------------------------------------------------
pdf("03_Output/Pairsplots_pairsfun.pdf", width = 15, height = 7)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]
  
  color_scheme_set("purple")
  # adding title to each plot 
  x <- pairs(stan_toutputs[[i]], pars = c("r", "k"), gap = 0, 
             main = paste0("Pairs plot_", strain), pch = 16, cex = 0.5) 
  
  print(x)
}

dev.off()


# mcmc_pairs fun ----------------------------------------------------------------


pdf("bacterial_comms_output/mcmcpairs_5_5_26_rk.pdf", width = 15, height = 7, onefile = TRUE)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]
  pst_cp2 <- as.array(stan_toutputs[[i]])
  
  color_scheme_set("brightblue")
  # adding title to each plot 
  x <- mcmc_pairs(pst_cp2, pars = c("r", "k"),
                  diag_fun = "hist", off_diag_fun = "hex",
                  condition = pairs_condition(chains = list(1, 2:3)), 
                  grid_args = list(top = paste("Pairs plot:", strain))) 
  
  print(x)
}

dev.off()

# r values bxpl -----------------------------------------------------------------

list_r <- lapply(names(stan_toutputs), function(nm) {
  # 
  
  draws <- as_draws_df(stan_toutputs[[nm]]) %>%
    select(r) %>% # select the r value in each spp 
    mutate(ID = nm) # keep the name 
  return(draws)
})

# Bind all the rows in the object and make two columns one for each variable 
df_boxplot <- bind_rows(list_r) %>% 
separate(ID, into = c("Strain", "Temperature"), sep = "_T") 

df_boxplot <- df_boxplot %>%
   
  # generate a new "sampling_times column 
  mutate(
      sampling_times = case_when(
      str_detect(Temperature, "_6t")  ~ 6, # for the CH29 special cases 
      str_detect(Temperature, "_10t") ~ 10,
      TRUE                           ~ 10
    ),
    
    # to clean the temperature column 
    Temperature = str_replace_all(Temperature, "_6t|_10t", "")
  )

df_boxplot$Temperature <- as.numeric(df_boxplot$Temperature)
bxplot_rdf <- as.data.frame(df_boxplot)

ggplot(bxplot_rdf, aes(y = r, 
                x = Strain, 
                fill = factor(Temperature))) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~ Strain, scale="free", labeller = as_labeller(c("CH111" = "Bacillus thuringiensis",  "CH149a" = "Micrococcus luteus",  
                                                              "CH154a" = "Staphylococcus shinii", "CH161d" = "Bacillus infantis",
                                                              "CH23" = "Bacillus altitudinis",   "CH29" = "Corunebacterium sp", 
                                                              "CH447" = "Priestia megaterium",  "CH450" = "Metabacillus indicus",
                                                              "CH90" = "Bacillus atrophaeus", "CH99b" = "Staphylococcus arlettae"))) + 
  scale_fill_manual(values = c("rosybrown1", "rosybrown3", "pink4"))+
  scale_y_continuous(name = "r values") +
  scale_x_discrete(name = "Strain") +
  ggtitle("r values by temperature") +
  theme(plot.title = element_text(hjust = 0.5))

# k values bxpl -----------------------------------------------------------------

list_k <- lapply(names(stan_toutputs), function(nm) {
  # 
  
  drawsk <- as_draws_df(stan_toutputs[[nm]]) %>%
    select(k) %>% # select the r value in each spp 
    mutate(ID = nm) # keep the name 
  return(drawsk)
})

# Bind all the rows in the object and make two columns one for each variable 
df_boxplotk <- bind_rows(list_k) %>% 
  separate(ID, into = c("Strain", "Temperature"), sep = "_T") 

df_boxplotk <- df_boxplotk %>%
  
  # generate a new "sampling_times column 
  mutate(
    sampling_times = case_when(
      str_detect(Temperature, "_6t")  ~ 6, # for the CH29 special cases 
      str_detect(Temperature, "_10t") ~ 10,
      TRUE                           ~ 10
    ),
    
    # to clean the temperature column 
    Temperature = str_replace_all(Temperature, "_6t|_10t", "")
  )

df_boxplotk$Temperature <- as.numeric(df_boxplotk$Temperature)
bxplot_kdf <- as.data.frame(df_boxplotk)

ggplot(bxplot_kdf, aes(y = k, 
                       x = Strain, 
                       fill = factor(Temperature))) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~ Strain, scale="free", labeller = as_labeller(c("CH111" = "Bacillus thuringiensis",  "CH149a" = "Micrococcus luteus",  
                                                              "CH154a" = "Staphylococcus shinii", "CH161d" = "Bacillus infantis",
                                                              "CH23" = "Bacillus altitudinis",   "CH29" = "Corunebacterium sp", 
                                                              "CH447" = "Priestia megaterium",  "CH450" = "Metabacillus indicus",
                                                              "CH90" = "Bacillus atrophaeus", "CH99b" = "Staphylococcus arlettae"))) + 
  scale_fill_manual(values = c("tomato3", "tomato4", "red4"))+
  scale_y_continuous(name = "r values") +
  scale_x_discrete(name = "Strain") +
  ggtitle("k values by temperature") +
  theme(plot.title = element_text(hjust = 0.5)) 

# stan.hist ----------------------------------------------------------------------

pdf("03_Output/bayesplot_rk_plots/Histogram_rk.pdf", width = 10, height = 5 , onefile = TRUE)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]
  
  color_scheme_set("brightblue")
  # adding title to each plot 
  x <- stan_hist(stan_toutputs[[i]], pars = c("r", "k", "sigma")) + 
    ggtitle(paste("Pairs plot:", strain))
                  
  
  print(x)
}

dev.off()

#### POSTERIOR DISTRIBUTION VISUALIZATION ####
# plot density --------------------------------------------------------------------
pdf("bacterial_comms_output/densityplot_6_5_26.pdf", width = 15, height = 7 , onefile = TRUE)

for (i in 1:length(stan_toutputs)){
  
  strain <- names(stan_toutputs)[i]
  posterior <- as.array(stan_toutputs[[i]])
  

    p_dens <- mcmc_areas(
    posterior,
    pars = c("r", "k"),
    prob = 0.95,       # 95% int
    prob_outer = 0.99, # 99% line
    point_est = "median"
  ) +
    ggtitle(paste("Posterior distribution:", names(stan_toutputs)[1])) +
    theme_minimal()
  
  print(p_dens)
}
dev.off()

# mcmc_dens -------------------------------------------------------------------

pdf("bacterial_comms_output/densityplot_6_5_26.pdf", width = 10, height = 5 , onefile = TRUE)

for (i in 1:length(stan_toutputs)){
  
  strain <- names(stan_toutputs)[i]
  f <- as_draws(stan_toutputs[[i]])
  
  color_scheme_set("mix-teal-pink")
  
  g <- mcmc_dens_overlay(f, pars = c("r", "k"),
                    facet_args = list(nrow = 2)) +
    facet_text(size = 14) + 
    ggtitle(paste("Density", names(stan_toutputs)[i])) +
    theme_minimal()
  
  print(g)
}
dev.off()


# mcmc_areas -----------------------------------------------------------------

pdf("bacterial_comms_output/mcmc_areas_6_5_26.pdf", width = 10, height = 5 , onefile = TRUE)

for (i in 1:length(stan_toutputs)){
  
  strain <- names(stan_toutputs)[i]
  f <- as_draws_df(stan_toutputs[[i]])
  
  color_scheme_set("mix-teal-pink")
  
  g <- mcmc_areas(f, regex_pars = c("r", "k"),  prob = 0.8) +
    labs(
      title = paste("Posterior distributions", strain),
      subtitle = "with medians and 80% intervals"
    )
  
  print(g)
}
dev.off()

