#rm(list=ls())
library(deSolve)
library(DAISIE)
library(pracma)
###############################################################################
### fonction to calculate the likelihood of an island with no survivants descendants
###############################################################################

### Using D-E approach



# Define system of equations for interval [t0, tp]
DAISIE_DE_logp0 <- function(datalist, pars1)
{
  t0 <- datalist[[1]]$island_age
  tp <- 0
  parameters <- pars1
  interval0 <- function(t, state, parameters) {
      with(as.list(c(state, parameters)),
           {
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[2]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD0, dDm, dE1))
           })
    }

    # Set initial conditions
      initial_conditions0 <- c(D0 = 1, Dm = 0, E1 = 0)

    # Time sequence for interval [t0, tp]
        time0 <- c(tp, t0)

    # Solve the system for interval [t0, tp]
    solution0 <- ode(y = initial_conditions0,
               times = time0,
               func = interval0,
               parms = parameters,
               method = "lsodes",
               rtol= 1e-12, atol= 1e-12)

    # Extract log-likelihood
      L0 <- solution0[, "D0"][[2]]
      logL0b <- log(L0)
      return(logL0b)
}
#DAISIE_DE_logp0(datalist, pars1)
