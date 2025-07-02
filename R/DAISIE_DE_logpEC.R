#' @name DAISIE_DE_logpEC
#' @title Function to calculate the likelihood of observing an endemic lineage
#' with fixed colonization time. This is valid for infinite K according to the
#' DE equations.
#' @description This function calculates the log-likelihood of observing an
#' endemic lineage with fixed colonization time. This is valid for infinite K
#' according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @examples
#'
#' library(DAISIE)
#' data("Galapagos_datalist")
#' datalist <- Galapagos_datalist
#'
#' i <- 5
#' brts <- datalist[[i]]$branching_times
#'
#' parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
#' brts <- datalist[[i]]$branching_times
#' missnumspec <- datalist[[i]]$missing_species
#' DAISIE_DE_logpEC(
#'   brts                    = brts,
#'   status                  = 2,
#'   parameter               = parameter,
#'   atol                    = 1e-15,
#'   missnumspec             = missnumspec,
#'   rtol                    = 1e-15,
#'   methode                 = "ode45"
#' )
#' @noRd

DAISIE_DE_logpEC <- function(brts,
                             status,
                             parameter,
                             missnumspec,
                             atol  = 1e-15,
                             rtol  = 1e-15,
                             methode                 = "ode45",
                             rcpp_methode = "odeint::bulirsch_stoer",
                             use_Rcpp = 0) {
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  ti <- sort(brts)
  ti <- ti[1:(length(ti)-2)]
  time2 <- c(t2, t1)
  time3 <- c(t1, t0)

  # Initial conditions
  number_of_species <- length(brts) - 1
  sampling_fraction <- number_of_species / (missnumspec + number_of_species)
  
  if (status == 2) {
  initial_conditions1 <- get_initial_conditions2(status = status,
                                                brts = brts,
                                                sampling_fraction = sampling_fraction
                                                )
  
  solution0 <- deSolve::ode(y = initial_conditions1,
                            times = c(0, ti),
                            func = interval2_EC,
                            parms = parameter,
                            method = methode,
                            rtol = rtol,
                            atol = atol)

  
  # Time sequences for interval [t2, tp]
  times <- rbind(c(0, ti[1:(length(ti) - 1)]), ti)
  
  for (idx in 1:length(ti)) {
    # Time sequence idx in interval [t2, tp]
    time1 <- times[, idx]
    
    # Solve the system for interval [t2, tp]
    solution1 <- solve_branch( interval_func = interval2_EC,
                                   initial_conditions = initial_conditions1,
                                   time = time1,
                                   parameter = parameter,
                                   methode = methode,
                                   rcpp_methode = rcpp_methode,
                                   atol = atol,
                                   rtol = rtol,
                                   use_Rcpp = use_Rcpp)
    
    
    

    
    # Initial conditions
    initial_conditions1 <- c(DE = parameter[1] * solution0[, "DE"][idx + 1] * solution1[, "DE"][2],
                             DM3 = 0, E = solution0[, "E"][idx + 1], DA3 = 1)
  }
  
  # Initial conditions
  initial_conditions2 <- c(DE = initial_conditions1["DE"][[1]],
                           DM2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti) + 1],
                           DM3 = solution0[, "DM3"][length(ti) + 1],
                           E = initial_conditions1["E"][[1]],
                           DA3 = solution0[, "DA3"][length(ti) + 1])
  
  
  # Solve the system for interval [t2, tp]
  solution2 <-  solve_branch(interval_func = interval2_ES,
                              initial_conditions = initial_conditions2,
                              time = time2,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp
                            )

  # Initial conditions
  initial_conditions3 <- get_initial_conditions4(
                                                  status = status,
                                                  solution = solution2,
                                                  parameter = parameter
                                                )
  
  
                solution3 <- solve_branch(
                                          interval_func = interval4,
                                          initial_conditions = initial_conditions3,
                                          time = time3,
                                          parameter = parameter,
                                          methode = methode,
                                          rcpp_methode = rcpp_methode,
                                          atol = atol,
                                          rtol = rtol,
                                          use_Rcpp = use_Rcpp)
                

     }else if (status == 6) {
    
    initial_conditions1 <- get_initial_conditions2(status = status,
                                                   brts = brts,
                                                   sampling_fraction = sampling_fraction
    )
    solution0 <- deSolve::ode(y = initial_conditions1,
                              times = c(0, ti),
                              func = interval2_EC,
                              parms = parameter,
                              method = methode,
                              rtol = rtol,
                              atol = atol)
    
    # Time sequences for interval [t2, tp]
    times <- rbind(c(0, ti[1:(length(ti) - 1)]), ti)
    
    for (idx in 1:length(ti)) {
      # Time sequence idx in interval [t2, tp]
      time1 <- times[, idx]
      
      # Solve the system for interval [t2, tp]
      solution1 <- deSolve::ode(y = initial_conditions1,
                                times = time1,
                                func = interval2_EC,
                                parms = parameter,
                                method = methode,
                                rtol = rtol,
                                atol = atol)
      
      # Initial conditions
      initial_conditions1 <- c(DE = parameter[1] * solution0[, "DE"][idx + 1] * solution1[, "DE"][2],
                               DM3 = 0, E = solution0[, "E"][idx + 1], DA3 = 1)
    }
    
    # Initial conditions
    initial_conditions2 <- c(DE = initial_conditions1["DE"][[1]],
                             DM1 = 0,
                             DM2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti)+1],
                             DM3 = solution0[, "DM3"][length(ti) + 1],
                             E = initial_conditions1["E"][[1]],
                             DA2 = 0,
                             DA3 = solution0[, "DA3"][length(ti) + 1])
    
    
    # Solve the system for interval [t2, tp]
    solution2 <- deSolve::ode(y = initial_conditions2,
                              times = time2,
                              func = interval3_ES,
                              parms = parameter,
                              method = methode,
                              rtol = rtol,
                              atol = atol)
    
    
  
    
    # Initial conditions
    initial_conditions3 <- get_initial_conditions4(
                                                  status = status,
                                                  solution = solution2,
                                                  parameter = parameter
                                                )
    
    # Solve the system for interval [t0, t1]
    solution3 <- solve_branch(
                              interval_func = interval4,
                              initial_conditions = initial_conditions3,
                              time = time3,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
     } else if (status == 3 & length(brts) > 2) {

       initial_conditions1 <- get_initial_conditions2(status = status,
                                                      brts = brts,
                                                      sampling_fraction = sampling_fraction
       )
       
       solution0 <- deSolve::ode(y = initial_conditions1,
                                 times = c(0, ti),
                                 func = interval2_EC,
                                 parms = parameter,
                                 method = methode,
                                 rtol = rtol,
                                 atol = atol)
       
       
       # Time sequences for interval [t2, tp]
       times <- rbind(c(0, ti[1:(length(ti) - 1)]), ti)
       
       for (idx in 1:length(ti)) {
         # Time sequence idx in interval [t2, tp]
         time1 <- times[, idx]
         
         # Solve the system for interval [t2, tp]
         solution1 <- solve_branch( interval_func = interval2_EC,
                                    initial_conditions = initial_conditions1,
                                    time = time1,
                                    parameter = parameter,
                                    methode = methode,
                                    rcpp_methode = rcpp_methode,
                                    atol = atol,
                                    rtol = rtol,
                                    use_Rcpp = use_Rcpp)
         
         
         
         
         
         # Initial conditions
         initial_conditions1 <- c(DE = parameter[1] * solution0[, "DE"][idx + 1] * solution1[, "DE"][2],
                                  DM3 = 0, E = solution0[, "E"][idx + 1], DA3 = 1)
       }
       
       # Initial conditions
       initial_conditions2 <- c(DE = initial_conditions1["DE"][[1]],
                                DM2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti) + 1],
                                DM3 = solution0[, "DM3"][length(ti) + 1],
                                E = initial_conditions1["E"][[1]],
                                DA3 = solution0[, "DA3"][length(ti) + 1])
       
       
       # Solve the system for interval [t2, tp]
       solution2 <-  solve_branch(interval_func = interval2_ES,
                                  initial_conditions = initial_conditions2,
                                  time = time2,
                                  parameter = parameter,
                                  methode = methode,
                                  rcpp_methode = rcpp_methode,
                                  atol = atol,
                                  rtol = rtol,
                                  use_Rcpp = use_Rcpp
       )
       
       # Initial conditions
       initial_conditions3 <- get_initial_conditions4(
         status = status,
         solution = solution2,
         parameter = parameter
       )
       
       
       solution3 <- solve_branch(
         interval_func = interval4,
         initial_conditions = initial_conditions3,
         time = time3,
         parameter = parameter,
         methode = methode,
         rcpp_methode = rcpp_methode,
         atol = atol,
         rtol = rtol,
         use_Rcpp = use_Rcpp)
       
       
     }
  
  # Extract log-likelihood
  Lk <- solution3[, "DA1"][[2]]
  logLkb <- log(Lk)
  return(logLkb)
}
