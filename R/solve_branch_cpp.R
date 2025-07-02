#' @keywords internal
solve_branch_cpp <- function(chosen_func,
                             initial_conditions,
                             time,
                             parameter,
                             trait_mainland_ancestor,
                             methode = "odeint::bulirsch_stoer",
                             atol = 1e-15,
                             rtol = 1e-15) {

  lambda_c <- parameter[[1]]
  mus      <- parameter[[2]]
  gammas   <- parameter[[4]]
  lambda_a <- parameter[[5]]

  solution <- cpp_solve(lambda_c,
                        lambda_a,
                        mus,
                        gammas,
                        chosen_func,
                        methode,
                        initial_conditions,
                        time,
                        atol,
                        rtol)

  res <- matrix(data = NA, nrow = 2, ncol = length(solution$states))
  res[1, ] <- initial_conditions
  res[2, ] <- solution$states

  return(res)
}
