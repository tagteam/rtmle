#' Function to simulate Competing Risk data.
#'
#' @title Simulate Competing risk Data
#'
#' @param N A double of the number of individuals
#' @param beta A 2X3 matrix with the effect of L0 and Treatment on two competing processes and censoring.
#' The columns represent Process 1, Process 2 and Censoring, while the rows represent L0 and Treatment.
#' @param eta Vector of  length 3 of shape parameters for the Weibull hazard with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of length 3 of scale parameters for the Weibull hazard.
#'
#' @return  Data frame containing the simulated competing risk data
#' @export
#'
#' @examples
#' sim_comp_risk_data(10)
sim_comp_risk_data <- function(N,
                               beta = NULL,
                               eta = rep(0.1,3),
                               nu = rep(1.1,3)
){

  at_risk <- function(x, k, L) {
  if(x == 0 | x == 1 | x == 2) return(1)
    else return(0)
  }
  if(is.null(beta)){
    beta <- matrix(0, ncol = 3, nrow = 2)
  }
  beta <- rbind(c(beta[1,],0), rep(0,4), c(beta[2,],0), rep(0,4))
  results <- sim_event_data(N, beta, c(eta,0), c(nu,0), at_risk, term_deltas = c(0,1,2))
  results <- results[, !c("L")]

  return(results)
}
