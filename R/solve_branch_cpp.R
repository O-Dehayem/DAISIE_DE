#' @keywords internal
solve_branch_cpp <- function(chosen_func,
                             initial_conditions,
                             time,
                             parameter,
                             methode = "odeint::bulirsch_stoer",
                             atol = 1e-15,
                             rtol = 1e-15) {

  lambda_c <- parameter[[1]]
  mu      <- parameter[[2]]
  gamma   <- parameter[[4]]
  lambda_a <- parameter[[5]]

  solution <- cpp_solve(lambda_c,
                        lambda_a,
                        mu,
                        gamma,
                        chosen_func,
                        methode,
                        initial_conditions,
                        time,
                        atol,
                        rtol)

  res <- matrix(data = NA, nrow = 2, ncol = length(solution$states))
  res[1, ] <- initial_conditions
  res[2, ] <- solution$states
  colnames(res) <- names(initial_conditions)
  return(res)
}
