#' Function to simulate data from an T2D diabetes setting. 3 different types of events are simulated, chosen to represent
#' Death (1), Censoring(2) and Change in Covariate Process(3). Death and Censoring are terminal events and Change in
#' Covariate Process can occur once.
#' The intensities of the various events depend upon previous events and the pre specified \eqn{nu} and \eqn{eta}
#' parameters. The dependence on previous events is controlled by parameters chosen so that a large baseline covariate
#' (L0) increases the probability of covariate change (L = 1). Initial treatment (A0 = 1) reduces the risk of a change
#' in the covariate process (L= 1) and the risk of death. A large baseline covariate (L0) increase the risk of death.
#' probability of operation (A = 1). Operation (A = 1) decreases the risk of death. The censoring process does not
#' depend on anything. The risk of death depend on change in the covariate process (L = 1) through beta_L_D.
#'
#'
#' @title Simulate Data in a T2D Diabetes Setting
#'
#' @param N A double of the number of individuals
#' @param eta Vector of length 4 of shape parameters for the Weibull intensity with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of scale parameters for the Weibull hazard. Default is set to 0.1 for all events.
#' @param followup A censoring time. By default set to infinity.
#' @param beta_L_D Specifies how change in covariate process affects risk of death. Is by default set to 1.
#' @param cens Specifies whether you are at risk of being censored
#'
#' @return  Data frame containing the simulated data. There is a column for ID, time of event (Time),
#' event type (Delta), baseline covariate (L0), additional covariate (L) and  Treatment Process (A).
#' @export
#'
#' @examples
#' sim_data_setting2(10)
sim_data_setting2 <- function(N, eta = rep(0.1,4), nu = rep(1.1,4), followup = Inf, beta_L_D = 1,
                              cens = 1){

  at_risk <- function(i, L, A) {
    return(c(
      # You are only at risk for an operation if you have not had an operation yet
      0,
      # If you have not died yet or been censored yet, you are at risk for dying or being censored
      1,
      cens,
      # You are only at risk for a change in the covariate process if you have not had a change yet
      as.numeric(L[i] == 0)))
  }

  Time <- ID <- A <- NULL
  beta <- matrix(ncol = 4, nrow = 4)

  # No A
  beta[2,] <- 0; beta[,1] <- 0
  # L0 increases the probability of L = 1
  beta[1,4] <- 1
  # A0 = 1 reduces the risk of L = 1
  beta[3,4] <- -1
  # L0 increases the risk of death
  beta[1,2] <- 1
  # A0 reduces the risk of death
  beta[3,2] <- -1
  # L = 1 does not affect the intensity of L (the event occurs only once)
  beta[4,4] <- 0

  # Censorering does not depend on anything
  beta[,3] <- 0

  # Effect of L = 1 on risk of death
  beta[4,2] <- beta_L_D

  data <- sim_event_data(N, beta = beta, eta = eta, nu = nu, at_risk = at_risk, max_cens = followup)
  data[, A := NULL]

  # Incorporating followup time
  if(followup < Inf){

    # We first find the most recent event before followup for individuals that have been censored
    # and that have experienced an event
    data_cens <- data[ID %in% data[Time > followup]$ID & Time < followup]
    data_cens <- data_cens[,.SD[which.max(Time)], by=ID]

    # Now we are missing the individuals that have not yet experienced an event at the time of followup
    data_cens2 <- unique(data[Time > followup][!data[Time > followup]$ID %in% data_cens$ID, ], by = "ID")
    data_cens2$L <- 0

    # We combine the two cases
    data_cens <- rbind(data_cens, data_cens2)

    # We modify the time and event
    data_cens$Time <- followup
    data_cens$Delta <- NA

    # We combine and sort data
    data <- rbind(data[Time <= followup], data_cens)
    setkey(data, ID)
  }

  return(data)
}
