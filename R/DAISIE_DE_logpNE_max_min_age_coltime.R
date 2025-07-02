#' @name DAISIE_DE_logpNE_max_min_age_coltime
#' @title Function to calculate the likelihood of observing a non-endemic lineage on the island
#' with minimum and maximum times of colonization. This valid for infinite K according to the DE equations.
#' @description This function calculates the log-likelihood of observing a non-endemic lineage on an island
#' for which the exact colonization time is unknown, but the maximum and minimum times of colonization are
#' known. This is valid for infinite K according to the DE equations
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @noRd

#' @export

DAISIE_DE_logpNE_max_min_age_coltime <- function(brts,
                                                 status,
                                                 parameter,
                                                 methode,
                                                 rtol,
                                                 atol,
                                                 rcpp_methode = "odeint::bulirsch_stoer",
                                                 use_Rcpp = FALSE) {
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  time1 <- c(tp, t2)
  time2 <- c(t2, t1)
  time3 <- c(t1, t0)
  
  # Initial conditions
  initial_conditions1 <- get_initial_conditions2(
                                    status = status,
                                    brts = brts,
                                    sampling_fraction = 1
                                  )

  solution1 <- solve_branch(interval_func = interval2_NE,
                 initial_conditions = initial_conditions1,
                 time = time1,
                 parameter = parameter,
                 methode = methode,
                 rcpp_methode = rcpp_methode,
                 atol = atol,
                 rtol = rtol,
                 use_Rcpp = use_Rcpp)

  
  # Initial conditions
  initial_conditions2 <- c(DM1 = 0, DM2 = solution1[, "DM2"][[2]], E = solution1[, "E"][[2]], DA2 = 0)


  
  # Solve the system for interval [t1, tp]
  solution2 <-   solve_branch(interval_func = interval3_NE,
                              initial_conditions = initial_conditions2,
                              time = time2,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
  

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
