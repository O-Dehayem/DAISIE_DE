setwd('/home2/p306709/DAISIE_DE/results')
## load packages
library(DAISIE)




args <- as.numeric(commandArgs(TRUE))
i <- args[1]


i


ML_parameters_data_DI <- c(2.546591, 2.678781, 2, 0.009326754, 1.008583 )

# Simulated dataset

load ("/home2/p306709/DAISIE_DE/replicates_sim.RData")

Galapagos_sims_G <- output_c$Galapagos_sims_G






results <- DAISIE_ML(
  datalist = Galapagos_sims_G[[i]],
  initparsopt = c(2.546591, 2.678781, 2.6787812, 0.009326754, 1.008583),
  idparsopt = 1:5,
  parsfix = NULL,
  idparsfix = NULL,
  ddmodel = 0,
  methode = "ode45",
  CS_version = list(model = 1, function_to_optimize = 'DAISIE_DE'),
  equal_extinction = TRUE
)
results




output <- list(results = results)


save(output, file = paste( "data_G", i,".RData"))
