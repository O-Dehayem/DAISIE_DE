#' @name DAISIE_DE_logpES_max_min_age_coltime
#' @title Function to compute the likelihood of observing an endemic singleton lineage
#' on the island given the minimum and maximum colonization ages, valid for infinite K
#'according to the DE equations.
#' @description This function calculates the log-likelihood of observing an endemic singleton lineage on an island
#' for which the exact colonization time is unknown, but the maximum and minimum ages of colonization are given.
#' This is valid for infinite K according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @noRd

DAISIE_DE_logpES_max_min_age_coltime <- function(brts,
                                                 status,
                                                 missnumspec,
                                                 parameter,
                                                 methode,
                                                 rtol,
                                                 atol,
                                                 rcpp_methode = "odeint::bulirsch_stoer",
                                                 use_Rcpp = 0) {
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  time1 <- c(tp, t2)
  time2 <- c(t2, t1)
  time3 <- c(t1, t0)
  parameters <- pars1

  
  
  # Initial conditions
  
  initial_conditions1 <- get_initial_conditions2(status = status,
                                                 brts = brts,
                                                 sampling_fraction = sampling_fraction)
  
  solution1 <- solve_branch(interval_func = interval2_ES,
                            initial_conditions = initial_conditions1,
                            time = time1,
                            parameter = parameter,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol)

  
  # Initial conditions
  
  initial_conditions2 <- c(DE = solution1[, "DE"][[2]],
                           DM1 = 0,
                           DM2 = solution1[, "DM2"][[2]],
                           DM3 =  solution1[, "DM3"][[2]],
                           E =  solution1[, "E"][[2]],
                           DA2 = 0,
                           DA3 = solution1[, "DA3"][[2]])
  
  
  solution2 <- solve_branch(interval_func = interval3_ES,
                            initial_conditions = initial_conditions2,
                            time = time2,
                            parameter = parameter,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol)

  
  
  
  initial_conditions3 <- get_initial_conditions4(
    status = status,
    solution = solution2,
    parameter = parameter
  )
  
  
  solution3 <- solve_branch(interval_func = interval4,
                            initial_conditions = initial_conditions3,
                            time = time3,
                            parameter = parameter,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol,
                            use_Rcpp = use_Rcpp)
  
  # Extract log-likelihood
  L1 <- solution3[, "DA1"][[2]]
  logL1b <- log(L1)
  return(logL1b)
  
}

DAISIE_DE_logpES_max_min_age_coltime (brts,
                                       status = 9,
                                       missnumspec,
                                       parameter,
                                       methode,
                                       rtol,
                                       atol)

