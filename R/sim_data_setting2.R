#' Function to simulate data from a T2D diabetes setting. 3 different types of events are simulated, chosen to represent
#' Censoring(0), Death (1) and Change in Covariate Process(3). Death and Censoring are terminal events and Change in
#' Covariate Process can occur once.
#' The intensities of the various events depend upon previous events and the pre specified \eqn{\nu} and \eqn{\eta}
#' parameters. The dependence on previous events is controlled by parameters chosen so that a large baseline covariate
#' (L0) increases the probability of covariate change (L = 1). Initial treatment (A0 = 1) reduces the risk of a change
#' in the covariate process (L= 1) and the risk of death. A large baseline covariate (L0) increase the risk of death.
#' Censoring does not depend on anything. The risk of death depends on change in the covariate process (L = 1) through beta_L_D.
#'
#'
#' @title Simulate Data in a T2D Diabetes Setting
#'
#' @param N A double of the number of individuals
#' @param eta Vector of length 4 of shape parameters for the Weibull intensity with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of scale parameters for the Weibull hazard. Default is set to 0.1 for all events.
#' @param followup A maximal censoring time. By default set to infinity.
#' @param beta_L_D Specifies how change in the covariate process affects risk of death. Is by default set to 1.
#' @param beta_L0_L Specifies how change in the covariate process affects risk of death. Is by default set to 1.
#' @param beta_A0_D Specifies how baseline treatment affects risk of death. Is by default set to 0.
#' @param beta_A0_L Specifies how baseline treatment affects risk of T2D. Is by default set to 0.
#' @param cens Specifies whether you are at risk of being censored
#' @param sex A TRUE/FALSE indicating whether there should be an additional binary covariate L1, representing e.g. sex.
#' If this the case, sex = TRUE, an additional row in the beta matrix specifies what effect the covariate has on the
#' intensities of the varies events. If this row is not specified the deafualt effect is 0.
#'
#' @return  Data frame containing the simulated data. There is a column for ID, time of event (Time),
#' event type (Delta), baseline covariate (L0) and additional covariate (L).
#' @export
#'
#' @examples
#' sim_data_setting2(10)
sim_data_setting2 <- function(N, eta = rep(0.1,4), nu = rep(1.1,4), followup = Inf,
                              beta_L0_L = 1, beta_L_D = 1, beta_A0_D = 0,
                              beta_A0_L = 0, cens = 1, sex = TRUE,
                              beta_L0_D = 1,){

  at_risk <- function(i, L, A) {
    return(c(
      cens,
      # If you have not died yet or been censored yet, you are at risk for dying or being censored
      1,
      # You are never at risk for an operation
      0,
      # You are only at risk for a change in the covariate process if you have not experienced a change yet
      as.numeric(L[i] == 0)))
  }

  Time <- ID <- A <- NULL
  beta <- matrix(ncol = 4, nrow = 5)

  # No A
  beta[4,] <- 0; beta[,3] <- 0
  # How L0 affects the probability of L = 1
  beta[1,4] <- beta_L0_L
  # How A0 = 1 affects the risk of L = 1
  beta[2,4] <- beta_A0_L
  # L0 increases the risk of death
  beta[1,2] <- beta_L0_D
  # How A0 affects the risk of death
  beta[2,2] <- beta_A0_D
  # L = 1 does not affect the intensity of L (the event occurs only once)
  beta[3,4] <- 0
  # Censorering does not depend on anything
  beta[,1] <- 0
  # Effect of L = 1 on risk of death
  beta[3,2] <- beta_L_D
  # Sex has no effect
  beta[5,] <- 0

  data <- sim_event_data(N, beta = beta, eta = eta, nu = nu, at_risk = at_risk,
                         max_cens = followup, sex = TRUE)
  data[, A := NULL]

  return(data)
}
