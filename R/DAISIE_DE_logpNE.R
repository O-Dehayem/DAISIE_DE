
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
#' i <- 3
#' brts <- datalist[[i]]$branching_times
#'
#' parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
#' brts <- datalist[[i]]$branching_times
#' DAISIE_DE_logpNE(
#'   brts                    = brts,
#'   status                  = 4,
#'   parameter               = parameter,
#'   atol                    = 1e-15,
#'   rtol                    = 1e-15,
#'   methode                 = "ode45"
#' )
DAISIE_DE_logpNE <- function(brts,
                             status,
                             parameter,
                             sampling_fraction = NA,
                             atol = 1e-15,
                             rtol = 1e-15,
                             methode = "ode45",
                             rcpp_methode = "odeint::bulirsch_stoer",
                             use_Rcpp = 0) {
  
  # Unpack times from brts
  t0 <- brts[1]
  tmax <- brts[2]
  t1 <- brts[2]
  tp <- 0
  
  # Define time intervals
  time2 <- c(tp, t1)
  time3 <- c(tp, tmax)
  time4 <- c(tmax, t0)
  
  # Initialize solution4 based on status
  if (status == 4) {
    initial_conditions2 <- get_initial_conditions2(
      status = status,
      brts = brts,
      sampling_fraction = sampling_fraction
    )
    
    solution2 <- solve_branch(
      interval_func = interval2_NE,
      initial_conditions = initial_conditions2,
      time = time2,
      parameter = parameter,
      methode = methode,
      rcpp_methode = rcpp_methode,
      atol = atol,
      rtol = rtol,
      use_Rcpp = use_Rcpp
    )
    
    initial_conditions4 <- get_initial_conditions4(
      status = status,
      solution = solution2,
      parameter = parameter
    )
    
    solution4 <- solve_branch(
      interval_func = interval4,
      initial_conditions = initial_conditions4,
      time = time4,
      parameter = parameter,
      methode = methode,
      rcpp_methode = rcpp_methode,
      atol = atol,
      rtol = rtol,
      use_Rcpp = use_Rcpp
    )
    
  } else if (status == 1) {
    
    initial_conditions3 <- get_initial_conditions3(
      status = status,
      sampling_fraction = sampling_fraction
    )
    
    solution3 <- solve_branch(
      interval_func = interval3_NE,
      initial_conditions = initial_conditions3,
      time = time3,
      parameter = parameter,
      methode = methode,
      rcpp_methode = rcpp_methode,
      atol = atol,
      rtol = rtol
    )
    
    initial_conditions4 <- get_initial_conditions4(
      status = status,
      solution = solution3,
      parameter = parameter
    )
    
    solution4 <- solve_branch(
      interval_func = interval4,
      initial_conditions = initial_conditions4,
      time = time4,
      parameter = parameter,
      methode = methode,
      rcpp_methode = rcpp_methode,
      atol = atol,
      rtol = rtol,
      use_Rcpp = use_Rcpp
    )
    
  } else {
    stop("Unsupported status. Currently, only status 1 and 4 are handled.")
  }
  
  # Extract log-likelihood
  if (!"DA1" %in% colnames(solution4)) {
    stop("Expected column 'DA1' not found in solution4.")
  }
  
  Lk <- solution4[, "DA1"][[2]]
  logLkb <- log(Lk)
  return(logLkb)
}
