test_that("logpNE_max_min_age_coltime", {
  
  if (requireNamespace("DAISIE", quietly = TRUE)) {
    
    data("Galapagos_datalist", package = "DAISIE")
    datalist <- Galapagos_datalist
    brts <- c(5, 4, 3)
    
    parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
    
    res1 <- DAISIE_DE_logpNE_max_min_age_coltime(
      brts = brts,
      status = 8,
      parameter = parameter,
      atol = 1e-15,
      rtol = 1e-15,
      methode = "ode45",
      rcpp_methode = "odeint::bulirsch_stoer",
      use_Rcpp = FALSE
    )
    
    res2 <- DAISIE:::DAISIE_loglik_CS_choice(
      pars1 = parameter,
      pars2 = c(100, 11, 0, 2),
      brts = brts,
      stac = 8,
      missnumspec = 0,
      datalist = datalist
    )
    

    
    testthat::expect_equal(res1, res2, tolerance = 1e-4)
  }
})
