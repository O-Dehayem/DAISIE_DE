
#' @keywords internal
interval2_NE <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    
    lambdac <- parameter[1]
    mu      <- parameter[2]
    gamma   <- parameter[4]
    lambdaa <- parameter[5]


    dDM2 <- -(lambdac + mu + gamma + lambdaa) * DM2

    dE <- mu - (mu + lambdac) * E + lambdac * E * E

    return(list(c(dDM2, dE)))
  })
}


interval2_ES <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    lambdac <- parameter[1]
    mu      <- parameter[2]
    gamma   <- parameter[4]
    lambdaa <- parameter[5]


    
    dDE <- -(lambdac + mu) * DE + 2 * lambdac * DE * E
    
    dDM2 <- -(lambdac + mu + gamma + lambdaa) * DM2 + (lambdaa * DE + 2 * lambdac * DE * E) * DA3
    
    dDM3 <- -(lambdac + mu + lambdaa) * DM3 + (mu + lambdaa * E + lambdac * E * E) * DA3
    
    dE <- mu - (mu + lambdac) * E + lambdac * E * E
    
    dDA3 <- -gamma * DA3 + gamma * DM3
    
    return(list(c(dDE, dDM2, dDM3, dE, dDA3)))
  })
}


interval2_EC <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    
    lambdac <- parameter[1]
    mu      <- parameter[2]
    gamma   <- parameter[4]
    lambdaa <- parameter[5]
    
 
    dDE <- -(lambdac + mu) * DE + 2 * lambdac * DE * E
    
    dDM3 <- -(lambdac + mu + lambdaa) * DM3 + (mu + lambdaa * E + lambdac * E * E) * DA3
    
    dE <- mu - (mu + lambdac) * E + lambdac * E * E
    
    dDA3 <- -gamma * DA3 + gamma * DM3
    
    return(list(c(dDE, dDM3, dE, dDA3)))
  })
}


#' @keywords internal
interval3_ES <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    
    lambdac <- parameter[1]
    mu      <- parameter[2]
    gamma   <- parameter[4]
    lambdaa <- parameter[5]

    
    dDE <- -(lambdac + mu) * DE + 2 * lambdac * DE * E 
    
    dDM1 <- -(lambdac + mu + lambdaa + gamma) * DM1 + gamma*DM2 + (mu + lambdaa * E + lambdac * E * E) * DA2
    
    dDM2 <- -(lambdac + mu + lambdaa) * DM2 + (mu + lambdaa * E + lambdac * E * E) * DA2 + (lambdaa * DE + 2 * lambdac * DE*E) * DA3
    
    dDM3 <- -(lambdac + mu + lambdaa) * DM3 + (mu + lambdaa * E + lambdac * E * E) * DA3
    
    dE <- mu - (mu + lambdac) * E + lambdac * E * E
    
    dDA2 <- -gamma * DA2 + gamma * DM2
    dDA3 <- -gamma * DA3 + gamma * DM3
    
    return(list(c(dDE, dDM1, dDM2, dDM3, dE, dDA2, dDA3)))
  })
}
#' @keywords internal
interval3_NE <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    
    lambdac <- parameter[1]
    mu      <- parameter[2]
    gamma   <- parameter[4]
    lambdaa <- parameter[5]
    

    
    dDM1 <- -(lambdac + mu + lambdaa + gamma) * DM1 + (mu + lambdaa * E + lambdac * E * E) * DA2 + gamma*DM2
    
    dDM2 <- -(lambdac + mu + lambdaa) * DM2 + (mu + lambdaa * E + lambdac * E * E) * DA2
    
    
    dE <- mu - (mu + lambdac) * E + lambdac * E * E
    
    dDA2 <- -gamma * DA2 + gamma * DM2

    return(list(c(dDM1, dDM2, dE, dDA2)))
  })
}

#' @keywords internal
interval4 <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    lambdac <- parameter[1]
    mu      <- parameter[2]
    gamma   <- parameter[4]
    lambdaa <- parameter[5]



    dDA1 <- -gamma * DA1 + gamma * DM1
 
    
    dDM1 <- -(lambdac + mu + lambdaa) * DM1 + (mu + lambdaa * E + lambdac * E * E) * DA1
    
    dE <- mu - (mu + lambdac) * E + lambdac * E * E
    
    return(list(c( dDA1, dDM1, dE)))
  })
}

