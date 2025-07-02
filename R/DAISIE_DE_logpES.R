
#' This function calculates the likelihood of observing a non-endemic lineage with specified species trait states,
#'
#'
#' @inheritParams default_params_doc
#'
#' @export
#'
#' @examples
#' library(DAISIE)
#' data("Galapagos_datalist")
#' datalist <- Galapagos_datalist
#'
#' i <- 6
#' brts <- datalist[[i]]$branching_times
#'
#' parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
#' brts <- datalist[[i]]$branching_times
#' missnumspec <- datalist[[i]]$missing_species
#' DAISIE_DE_logpES(
#'   brts                    = brts,
#'   status                  = 2,
#'   parameter               = parameter,
#'   atol                    = 1e-15,
#'   missnumspec             = missnumspec,
#'   rtol                    = 1e-15,
#'   methode                 = "ode45"
#' )
DAISIE_DE_logpES <- function(brts,
                             status,
                             parameter,
                             missnumspec,
                             atol  = 1e-15,
                             rtol  = 1e-15,
                             methode                 = "ode45",
                             rcpp_methode = "odeint::bulirsch_stoer",
                             use_Rcpp = FALSE) {
  
  
  
  # Unpack times from brts
  t0   <- brts[1]
  tmax <- brts[2]
  t1   <- brts[2]
  
  tp   <- 0
  
  # Time intervals
  
  time2 <- c(tp, t1)
  time3 <- c(tp, tmax)
  time4 <- c(tmax, t0)
  
  # Solve for interval [tp, t2] (stem phase)
  
  # Initial conditions
  
  number_of_species <- length(brts) - 1
  sampling_fraction <- number_of_species / (missnumspec + number_of_species)
  
  # Run appropriate sequence of intervals
  if (status == 2 & length(brts == 2)) {
    initial_conditions2 <- get_initial_conditions2(status = status,
                                                   brts = brts,
                                                   sampling_fraction = sampling_fraction)
    
    solution2 <- solve_branch(interval_func = interval2_ES,
                              initial_conditions = initial_conditions2,
                              time = time2,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
    
    initial_conditions4 <- get_initial_conditions4(status = status,
                                                   solution = solution2,
                                                   parameter = parameter)
    solution4 <- solve_branch(interval_func = interval4,
                              initial_conditions = initial_conditions4,
                              time = time4,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
    
  } else if (status == 5) {
    initial_conditions3 <- get_initial_conditions3(status = status, sampling_fraction = sampling_fraction)
    
    solution3 <- solve_branch(interval_func = interval3_ES,
                              initial_conditions = initial_conditions3,
                              time = time3,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
    
    initial_conditions4 <- get_initial_conditions4(status = status, solution = solution3, parameter = parameter)
    solution4 <- solve_branch(interval_func = interval4,
                              initial_conditions = initial_conditions4,
                              time = time4,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
  } else if (status == 3 & length(brts == 2)) {
    initial_conditions2 <- get_initial_conditions2(status = status,
                                                   brts = brts,
                                                   sampling_fraction = sampling_fraction)
    
    solution2 <- solve_branch(interval_func = interval2_ES,
                              initial_conditions = initial_conditions2,
                              time = time2,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
    
    initial_conditions4 <- get_initial_conditions4(status = status,
                                                   solution = solution2,
                                                   parameter = parameter)
    solution4 <- solve_branch(interval_func = interval4,
                              initial_conditions = initial_conditions4,
                              time = time4,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
  }
  
  # Extract log-likelihood from final solution
  Lk <- solution4[, "DA1"][[2]]
  logLkb <- log(Lk)
  return(logLkb)
}
