#' Function to simulate data from a Survival setting. The function simulates data corresponding to $N$ individuals that
#' are at risk for censoring (0) and an event (1).
#'
#' @title Simulate Survival Data
#'
#' @param N A double of the number of individuals
#' @param beta A 2X2 matrix with the effect of L0 and A0 on Censoring and Event.
#' The columns represent Censoring and Death, while the rows represent L0 and A0.
#' @param eta Vector of  length 2 of shape parameters for the Weibull hazard with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of length 2 of scale parameters for the Weibull hazard.
#'
#' @return  Data frame containing the simulated survival data
#' @export
#'
#' @examples
#' sim_surv_data(10)
sim_surv_data <- function(N,
                         beta = NULL,
                         eta = rep(0.1,2),
                         nu = rep(1.1,2)
){

  at_risk <- function(i, L, A) c(1,1,0,0)

  if(is.null(beta)){
    beta <- matrix(0, ncol = 2, nrow = 2)
  }
  beta <- rbind(c(beta[1,],0, 0), c(beta[2,],0, 0),  rep(0,4), rep(0,4))
  results <- sim_event_data(N, beta, c(eta,0,0), c(nu,0,0), at_risk)
  results <- results[, !c("L", "A")]

  return(results)
}
