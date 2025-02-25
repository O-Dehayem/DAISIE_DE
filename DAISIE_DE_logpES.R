#rm(list=ls())
library(deSolve)
library(DAISIE)
library(pracma)
###############################################################################
### fonction to calculate the likelihood of observing an endemic singleton lineage at time t1
###############################################################################

### Using D-E approach

# pars2[1] corresponds to lx = length of ODE variable x
# pars2[2] = 11: linear dependence in speciation rate and in immigration rate
# pars2[3] = 0: corresponds to conditioning on island age
# pars2[4] = 1: sets that parameters and likelihood should be printed

# pars1[1] corresponds to the Cladogenesis rate
# pars1[2] corresponds to the Extinction rate of endemic lineages
# pars1[3] corresponds to the Extinction rate of non-endemic lineages
# pars1[4] = corresponds to the Colonization rate
# pars1[5] = corresponds to the Anagenesis rate


Likelihood_ES_lineage <- function(datalist, i, pars1) {
  t0 <- datalist[[1]]$island_age
  t1 <- datalist[[i]]$branching_times[2]
  tp <- 0
  parameters <- pars1
  # Define system of equations for interval [t1, tp]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dD1 <- -(pars1[1] + pars1[2]) * D1 + 2 * pars1[1] * D1 * E1
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0 
      dDM <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * DM + (pars1[5] * D1 + 2 * pars1[1] * D1 * E1) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD1, dD0, dDm, dDM, dE1))
    })
  }
  
  # Define system of equations for interval [t0, t1]
  interval2 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD0, dDm, dE1))
    })
  }
  
  # Initial conditions
  initial_conditions1 <- c(D1 = 1, D0 = 1, Dm = 0, DM = 0, E1 = 0)
  
  # Time sequence for interval [t1, tp]
  time1 <- c(tp, t1)
  
  # Solve the system for interval [t1, tp]
  solution1 <- ode(y = initial_conditions1,
                   times = time1,
                   func = interval1,
                   parms = parameters,
                   method = "lsodes",
                   rtol = 1e-12, atol = 1e-12)
  
  # Initial conditions
  initial_conditions2 <- c(D0 = pars1[4] * solution1[, "DM"][[2]],
                           Dm = pars1[4] * solution1[, "DM"][[2]],
                           E1 = solution1[, "E1"][[2]])
  
  # Time sequence for interval [t0, t1]
  time2 <- c(t1, t0)
  
  # Solve the system for interval [t0, t1]
  solution2 <- ode(y = initial_conditions2,
                   times = time2,
                   func = interval2,
                   parms = parameters, 
                   method = "lsodes",
                   rtol = 1e-12, atol = 1e-12)
  
  # Extract log-likelihood
  L1 <- solution2[, "D0"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}


