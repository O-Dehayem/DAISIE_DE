DAISIE_DE_loglik_CS <- function( parameter,
                                 pars2,
                                 datalist,
                                 methode = "lsodes",
                                 abstolint = 1e-15,
                                 reltolint = 1e-15)
{

  
  cond <- pars2[3]
  island_age <- datalist[[1]]$island_age
  
  if (length(parameter) == 5) {
    logp0 <- DAISIE_DE_logp0(datalist,
                             parameter,
                             atol = 1e-15,
                             rtol = 1e-15,
                             methode = "ode45",
                             rcpp_methode = "odeint::bulirsch_stoer",
                             use_Rcpp = 0)
    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
      numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik <- datalist[[1]]$not_present * logp0
      numimm <- datalist[[1]]$not_present + length(datalist) - 1
    }
    logcond <- (cond == 1) * log(1 - exp(numimm * logp0))
    for (i in 2:length(datalist)) {
      datalist[[i]]$type1or2 <- 1
    }
  } else {
    numimm <- datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 + length(datalist) - 1
    numimm_type2 <- length(which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2))
    numimm_type1 <- length(datalist) - 1 - numimm_type2
    if (!is.na(parameter[11])) {
      if (parameter[11] < numimm_type2 / numimm | parameter[11] > (1 - numimm_type1 / numimm)) {
        return(-Inf)
      }
      datalist[[1]]$not_present_type2 <- max(0, round(parameter[11] * numimm) - numimm_type2)
      datalist[[1]]$not_present_type1 <- numimm - (length(datalist) - 1) - datalist[[1]]$not_present_type2
    }
    logp0_type1 <- DAISIE_DE_logp0(datalist,
                                   parameter,
                                   atol = 1e-15,
                                   rtol = 1e-15,
                                   methode = "ode45",
                                   rcpp_methode = "odeint::bulirsch_stoer",
                                   use_Rcpp = 0)
    logp0_type2 <- DAISIE_DE_logp0(datalist,
                                   parameter,
                                   atol = 1e-15,
                                   rtol = 1e-15,
                                   methode = "ode45",
                                   rcpp_methode = "odeint::bulirsch_stoer",
                                   use_Rcpp = 0)
    loglik <- datalist[[1]]$not_present_type1 * logp0_type1 + datalist[[1]]$not_present_type2 * logp0_type2
    logcond <- (cond == 1) * log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) * logp0_type1 +
                                           (datalist[[1]]$not_present_type2 + numimm_type2) * logp0_type2))
  }
  
  loglik <- loglik - logcond
  vec_loglikelihood <- c()
  
  for (i in 2:length(datalist)) {
    status <- datalist[[i]]$stac
    brts <- datalist[[i]]$branching_times
    missnumspec <- datalist[[i]]$missing_species
    
    if (status == 1 ||status == 4) {
      loglikelihood <- DAISIE_DE_logpNE(brts,
                                        status,
                                        parameter,
                                        missnumspec,
                                        atol  = 1e-15,
                                        rtol  = 1e-15,
                                        methode                 = "ode45",
                                        rcpp_methode = "odeint::bulirsch_stoer",
                                        use_Rcpp = 0)
    } else if (status == 2 && length(brts) == 2 || status == 3 && length(brts) == 2 || status == 5 && length(brts) == 2 || status == 6) {
  
        loglikelihood <- DAISIE_DE_logpES(brts,
                                          status,
                                          parameter,
                                          missnumspec,
                                          atol  = 1e-15,
                                          rtol  = 1e-15,
                                          methode                 = "ode45",
                                          rcpp_methode = "odeint::bulirsch_stoer",
                                          use_Rcpp = 0)
    } else if (status == 2 && length(brts) > 2 || status == 3 && length(brts) > 2 || status == 6) {
      
      loglikelihood <- DAISIE_DE_logpEC(brts,
                                        status,
                                        parameter,
                                        missnumspec,
                                        atol  = 1e-15,
                                        rtol  = 1e-15,
                                        methode                 = "ode45",
                                        rcpp_methode = "odeint::bulirsch_stoer",
                                        use_Rcpp = 0)
    }  else if (status == 8) {
      loglikelihood <- DAISIE_DE_logpNE_max_min_age_coltime(brts,
                                                            status,
                                                            parameter,
                                                            missnumspec,
                                                            atol  = 1e-15,
                                                            rtol  = 1e-15,
                                                            methode                 = "ode45",
                                                            rcpp_methode = "odeint::bulirsch_stoer",
                                                            use_Rcpp = 0)
    } else if (status == 9) {
      loglikelihood <- DAISIE_DE_logpES_max_min_age_coltime(brts,
                                                            status,
                                                            parameter,
                                                            missnumspec,
                                                            atol  = 1e-15,
                                                            rtol  = 1e-15,
                                                            methode                 = "ode45",
                                                            rcpp_methode = "odeint::bulirsch_stoer",
                                                            use_Rcpp = 0)
    } else {
      stop("Unknown status value: ", status)
    }
    
    vec_loglikelihood <- c(vec_loglikelihood, loglikelihood)
    
    DAISIE:::print_parameters_and_loglik(
      pars = c(status, parameter),
      loglik = loglikelihood,
      verbose = pars2[4],
      parnames = c("lambda^c", "mu", "K", "gamma", "lambda^a", "prob_init_pres"),
      type = 'clade_loglik'
    )
  }
  loglik <- sum(vec_loglikelihood) + loglik
  return(loglik)
}

