#' Function to simulate survival data.
#'
#' @title Simulate Survival Data
#'
#' @param N A double of the number of individuals
#' @param beta A 2X2 matrix with the effect of L0 and Treatment on Death and Censoring.
#' The columns represent Death and Censoring, while the rows represetn L0 and Treatment.
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

  at_risk <- function(i) c(0,1,1,0)

  if(is.null(beta)){
    beta <- matrix(0, ncol = 2, nrow = 2)
  }
  beta <- rbind(c(0, beta[1,],0), rep(0,4), c(0, beta[2,],0), rep(0,4))
  results <- vec_sim_event_data(N, beta, c(0,eta,0), c(0,nu,0), at_risk)
  results <- results[, !c("L")]

  return(results)
}
