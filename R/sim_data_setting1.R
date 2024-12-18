#' Function to simulate data from an Operation setting. 4 different types of events are simulated, chosen to represent
#' Operation (0),Death (1), Censoring(2) and Change in Covariate Process(3). Death and Censoring are terminal events
#' and Operation and Change in Covariate Process can occur once.
#' The intensities of the various events depend upon previous events and the pre specified \eqn{nu} and \eqn{eta}
#' parameters. The dependence on previous events is controlled by parameters chosen so that a large baseline covariate
#' (L0) increases the probability of operation(A = 1) and covariate change (L = 1). A change in the covariate process
#' (L= 1) increases the probability of operation (A = 1). Operation (A = 1)decreases the risk of death. A large baseline
#' covariate (L0), and a change in the covariate process(L = 1) increase the risk of death. The censoring process does
#' not depend on anything.
#'
#' @title Simulate Operation Data in Setting 1
#'
#' @param N A double of the number of individuals
#' @param eta Vector of length 4 of shape parameters for the Weibull intensity with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of scale parameters for the Weibull hazard. Default is set to 0.1 for all events.
#' @param followup A censoring time. By default set to infinity.
#'
#' @return  Data frame containing the simulated data. There is a column for ID, time of event (Time),
#' event type (Delta), baseline covariate (L0), additional covariate (L) and  Treatment Process (A).
#' @export
#'
#' @examples
#' sim_data_setting1(10)
sim_data_setting1 <- function(N, eta = rep(0.1,4), nu = rep(1.1,4), followup = Inf){

  Time <- A0 <- ID <- NULL

  beta <- matrix(ncol = 4, nrow = 4)

  # A0 is 0
  beta[3,] <- 0
  # L0 increases the probability of L = 1 and A = 1
  beta[1,c(1,4)] <- c(1,1)
  # L=1 increases the probability of A = 1
  beta[4,1] <- 1
  # A=1 decreases the risk of death
  beta[2,2] <- -1
  # L0,L=1 increases the risk of death
  beta[c(1,4),2] <- 1
  # Censorering does not depend on anything
  beta[,3] <- 0

  # A = 1 does not affect the intensity of A (the event occurs only once)
  beta[2,1] <- 0
  # L = 1 does not affect the intensity of L (the event occurs only once)
  beta[4,4] <- 0
  # My assumption: A = 1 decreases the risk of L = 1
  beta[2,4] <- -0.5

  data <- sim_event_data(N, beta = beta, eta = eta, nu = nu)
  data[, A0 := NULL]

  # Incorporating followup time
  if(followup < Inf){

    # Identifying censored individuals
    cens <- unique(data[Time > followup], by = "ID")

    # Creating censored rows
    n <- nrow(cens)
    L0 <- cens$L0
    ID_cens <- cens$ID
    Time <- rep(followup, n)
    Delta <- rep("C", n)

    # Creating new rows for censored individuals
    cens_data <- data[ID %in% ID_cens & Time <= followup]

    # The A/L value right before censoring
    A <- numeric(n)
    A[ID_cens %in% cens_data$ID] <- cens_data[, max(A), by=ID]$V1
    A[! ID_cens %in% cens_data$ID] <- 0
    L <- numeric(n)
    L[ID_cens %in% cens_data$ID] <- cens_data[, max(L), by=ID]$V1
    L[! ID_cens %in% cens_data$ID] <- 0

    data <- rbind(data[Time <= followup], data.table(ID = ID_cens, Time, Delta, L0, L, A))
    setkey(data, ID)
  }

  return(data)
}
