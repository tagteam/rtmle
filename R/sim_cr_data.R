#' Function to simulate data from a Competing Risk setting. The function simulates data corresponding to $N$ individuals that
#' are at risk for mutually exclusive types of failure. 3 events can take place, one of which can be interpreted as censoring.
#'
#' @title Simulate Competing risk Data
#'
#' @param N A double of the number of individuals
#' @param beta A 2X3 matrix with the effect of L0 and A0 on the three processes.
#' The columns represent Process 1, Process 2 and Process 3, while the rows represent L0 and A0.
#' @param eta Vector of  length 3 of shape parameters for the Weibull hazard with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of length 3 of scale parameters for the Weibull hazard.
#'
#' @return  Data frame containing the simulated competing risk data. There is a column for ID, time of event (Time),
#' event type (Delta), baseline covariate (L0) and baseline Treatment (A0).
#' @export
#'
#' @examples
#' sim_cr_data(10)
sim_cr_data <- function(N,
                        beta = NULL,
                        eta = rep(0.1,3),
                        nu = rep(1.1,3)
                        ){

  at_risk <- function(i, L, A) c(1,1,1,0)

  if(is.null(beta)){
    beta <- matrix(0, ncol = 3, nrow = 2)
  }

  beta <- rbind(c(beta[1,],0), c(beta[2,],0), rep(0,4), rep(0,4))
  results <- sim_event_data(N, beta, c(eta,0), c(nu,0), at_risk, term_deltas = c(0,1,2))
  results <- results[, !c("L", "A")]

  return(results)
}
