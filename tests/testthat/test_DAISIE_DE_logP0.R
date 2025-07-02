test_that("DAISIE_logp0 is correct", {

  
  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist
  
  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
 
      
      res1 <- DAISIE:::DAISIE_DE_logp0(
        island_age = datalist[[1]]$island_age,
        pars1 = c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583),
        methode = "lsodes",
        reltolint = 1e-10,
        abstolint = 1e-10
      )
      
      res2 <- DAISIE_DE_logp0(datalist,
                              parameter,
                              atol = 1e-15,
                              rtol = 1e-15,
                              methode = "ode45",
                              rcpp_methode = "odeint::bulirsch_stoer",
                              use_Rcpp = 0
      )
      
      testthat::expect_equal(res1, res2, tolerance = 1e-6)
      
      res2 <- DAISIE_DE_logp0(datalist,
                              parameter,
                              atol = 1e-15,
                              rtol = 1e-15,
                              methode = "ode45",
                              rcpp_methode = "odeint::bulirsch_stoer",
                              use_Rcpp = 0
      )
      
      testthat::expect_equal(res1, res2, tolerance = 1e-6)
      

})
