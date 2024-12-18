#' Function to simulate data from an T2D diabetes setting. 3 different types of events are simulated, chosen to represent
#' Death (1), Censoring(2) and Change in Covariate Process(3). Death and Censoring are terminal events and Change in
#' Covariate Process can occur once.
#' The intensities of the various events depend upon previous events and the pre specified \eqn{nu} and \eqn{eta}
#' parameters. The dependence on previous events is controlled by parameters chosen so that a large baseline covariate
#' (L0) increases the probability of covariate change (L = 1). Initial treatment (A0 = 1) reduces the risk of a change
#' in the covariate process (L= 1) and the risk of death. A large baseline covariate (L0) increase the risk of death.
#' probability of operation (A = 1). Operation (A = 1) decreases the risk of death. The censoring process does not
#' depend on anything. Further two settings can be specified: "a") where the change in the covariate process (L = 1)
#' increases risk of death, and "b") where the risk of death does not depend on the change in the covariate process (L = 1).
#'
#'
#' @title Simulate Data in a T2D Diabetes Setting
#'
#' @param N A double of the number of individuals
#' @param eta Vector of length 4 of shape parameters for the Weibull intensity with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of scale parameters for the Weibull hazard. Default is set to 0.1 for all events.
#' @param followup A censoring time. By default set to infinity.
#' @param setting Specifies how change in covariate process affects risk of death. Can be set to either "a" or "b".
#' Is by default set to "a".
#'
#' @return  Data frame containing the simulated data. There is a column for ID, time of event (Time),
#' event type (Delta), baseline covariate (L0), additional covariate (L) and  Treatment Process (A).
#' @export
#'
#' @examples
#' sim_data_setting2(10)
sim_data_setting2 <- function(N, eta = rep(0.1,4), nu = rep(1.1,4), followup = Inf, setting = "a"){

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

  # 2a: L = 1 increases risk of death
  if(setting == "a") beta[4,2] <- 1

  # 2b: Risk of death does not depend on L
  else beta[4,2] <- 0

  data <- sim_event_data(N, beta = beta, eta = eta, nu = nu)
  data[, A := NULL]

  # Incorporating followup time
  if(followup < Inf){

    # Identifying censored individuals
    cens <- unique(data[Time > followup], by = "ID")

    n <- nrow(cens)
    A0 <- cens$A0
    L0 <- cens$L0
    ID_cens <- cens$ID
    Time <- rep(followup, n)
    Delta <- rep("C", n)

    # Creating new rows for censored individuals
    cens_data <- data[ID %in% ID_cens & Time <= followup]
    # The L value right before censoring
    L <- numeric(n)
    L[ID_cens %in% cens_data$ID] <- cens_data[, max(L), by = ID]$V1
    L[! ID_cens %in% cens_data$ID] <- 0

    data <- rbind(data[Time <= followup], data.table(ID = ID_cens, Time, Delta, L0, L, A0))
    setkey(data, ID)
  }

  return(data)
}
