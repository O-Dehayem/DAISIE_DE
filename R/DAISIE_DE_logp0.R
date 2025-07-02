#' testing fuction, for comparison with DAISIE
#' @description
#' This function compute the likelihood that all species that colonize the island
#' have gone extinct prior to the present.
#' @export
#' @inheritParams default_params_doc
#' @examples
#' #Load DAISIE package and data
#' library(DAISIE)
#' data("Galapagos_datalist")
#' datalist <- Galapagos_datalist
#'
#' parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
#' # Compute log‐likelihood under the DE‐trait model
#' DAISIE_DE_logp0(
#'   datalist            = datalist,
#'   parameter           = parameter,
#'   atol                = 1e-15,
#'   rtol                = 1e-15,
#'   methode             = "ode45"
#' )
#' @rawNamespace useDynLib(DAISIEde, .registration = TRUE)
#' @rawNamespace import(Rcpp)
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)

DAISIE_DE_logp0 <- function(datalist,
                            parameter,
                            atol = 1e-15,
                            rtol = 1e-15,
                            methode = "ode45",
                            rcpp_methode = "odeint::bulirsch_stoer",
                            use_Rcpp = FALSE) {
  

  t0 <- datalist[[1]]$island_age
  tp <- 0
  #########interval4 [t_p, t_0]
  
  initial_conditions40 <- c(DA1 = 1, DM1 = 0, E = 0)
  
  # Time sequence for interval [tp, t0]
  time4 <- c(tp, t0)
  
  # Solve the system for interval [tp, t1]
  solution4 <- solve_branch(interval_func = interval4,
                            initial_conditions = initial_conditions40,
                            time = time4,
                            parameter = parameter,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol,
                            use_Rcpp = use_Rcpp)
  
  # Extract log-likelihood
  Lk <- solution4[, "DA1"][[2]]
  logLkb <- log(Lk)
  return(logLkb)
}
