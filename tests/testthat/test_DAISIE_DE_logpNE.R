test_that("logpNE", {
  
  if (requireNamespace("DAISIE")) {
    data("Galapagos_datalist", package = "DAISIE")
    datalist <- Galapagos_datalist
    
    i <- 2
    brts <- datalist[[i]]$branching_times
   
    
    parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
    
    res1 <- DAISIE_DE_logpNE(brts,
                            status = 1,
                            parameter,
                            missnumspec,
                            atol  = 1e-15,
                            rtol  = 1e-15,
                            methode                 = "ode45",
                            rcpp_methode = "odeint::bulirsch_stoer",
                            use_Rcpp = 0)
    
    res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = c(2.546591, 2.678781, Inf, 0.009326754, 1.008583),
                                             pars2 = c(100, 11, 0, 2),
                                             brts = brts,
                                             stac = 1,
                                             missnumspec = 0,
                                             datalist = datalist)
    
    testthat::expect_equal(res1, res2)

  }
})

