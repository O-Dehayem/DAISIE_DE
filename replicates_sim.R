setwd('/home2/p306709/DAISIE_DE')
## load packages
library(DAISIE)
library(ape)
set.seed(1)



#Parameters use to simulate G dataset

ML_parameters_data_DI <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583 )
ML_parameters_data_DI[c(1,2,4,5)]

# Simulate replicates of island (G dataset)

data("Galapagos_datalist")

island_age <- Galapagos_datalist[[1]]$island_age
col_max <- island_age - 1e-05
replicates_to_run <- 14000


Galapagos_sims_G <- DAISIE_sim(
  time = island_age,
  M = 1000,
  pars = ML_parameters_data_DI,
  replicates = replicates_to_run,
  plot_sims = FALSE)


 




output_c <- list(Galapagos_sims_G = Galapagos_sims_G)


save(output_c, file = "replicates_sim.RData")
