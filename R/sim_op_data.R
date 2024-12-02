#' Function to simulate operation data.
#'
#' @title Simulate Operation Data
#'
#' @param N A double of the number of individuals
#' @param beta A 4X4 matrix of doubles for the effects on the intensities. The columns represent
#' the events Operation, Death, Censoring, and Covariate Change. The rows represent the baseline covariate \eqn{L0},
#' the event number \eqn{k}, treatment \eqn{A}, and the indicator \eqn{L} for change in the covariate process.
#' Default is set to 0.
#' @param eta Vector of length 4 of shape parameters for the Weibull hazard with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of scale parameters for the Weibull hazard. Default is set to 0.1 for all events.
#'
#' @return  Data frame containing the simulated operation data
#' @export
#'
#' @examples
#' sim_op_data(10)
sim_op_data <- function(N,
                        beta = NULL,
                        eta = rep(0.1,4),
                        nu = rep(1.1,4)){
  vec_sim_event_data(N, beta = beta, eta = eta, nu = nu)
}
