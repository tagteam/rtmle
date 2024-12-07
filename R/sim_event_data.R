#' Function to simulate event data. The number of events simulated corresponds
#' to the length of the eta and nu vector, and the number of columns in the beta matrix.
#' By default 4 different types of events are simulated: Operation (0), Death (1), Censoring(2)
#' and Change in Covariate Process(4). Death and Censoring are terminal events and Operation
#' and Change in Covariate Process can occur once.
#'
#' Different settings can be specified: Survival setting, competing risk setting and operation setting.
#' The different settings can be specified by the at_risk function. Default is the operation setting.
#'
#' @title Simulate Event Data
#'
#' @param N A double for the number of simulated individuals
#' @param beta A matrix of doubles for the effects on the intensities. The columns represent
#' the events Operation, Death, Censoring, and Covariate Change. The rows represent the baseline covariate \eqn{L0},
#' the event number \eqn{k}, treatment \eqn{A}, and the indicator \eqn{L} for change in the covariate process.
#' Default is set to 0.
#' @param eta Vector of shape parameters for the Weibull hazard with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of scale parameters for the Weibull hazard. Default is set to 0.1 for all events.
#' @param at_risk At risk function. Default is set to recurrent event setting. A survival
#' or competing risk setting can be specified as well. The \code{at_risk} function is an indicator of whether an individual is
#' at risk for a specific event. The function takes as input the event \eqn{x}, \eqn{k} and \eqn{L}.
#' @param term_deltas Terminal events. Default is set so that event 1 and 2 are terminal events.
#'
#' @return Data frame containing the simulated data. There is a column for ID, time of event (Time),
#' event type (Delta), baseline covariate (L0), additional covariate (L) and Treatment (A)
#' @export
#'
#' @examples
#' sim_event_data(N = 10)

sim_event_data <- function(N,               # Number of individuals
                    beta = NULL,            # Effects
                    eta = rep(0.1,4),       # Shape parameters
                    nu = rep(1.1,4),        # Scale parameters
                    at_risk = NULL,         # Function defining the setting
                    term_deltas = c(1,2)    # Terminal events
){

  if(is.null(beta)){
    beta <- matrix(0, nrow = 4, ncol = 4)
  }


  if(is.null(at_risk)){
    at_risk <- function(x, k, L) {
      # If you have not died yet or been censored yet, you are at risk for dying or being censored
      if(x == 1 | x == 2) return(1)
      # You are only at risk for an operation if you have not had an operation yet
      else if(x == 0) return(as.numeric(k == 0 | (k == 1 & L == 1)))
      # You are only at risk for a change in the covariate process if you have not had a change yet
      else return(as.numeric(L == 0))
      }
    }

  # Events
  x <- 1:ncol(beta)

  # Intensities
  phi_x <- function(x, k, L) {
    exp(L0 * beta[1,x] + k * beta[2,x] + A * beta[3,x] + L * beta[4,x])
  }

  lambda_x <- function(x, t, k, L) {
    at_risk(x - 1, k, L) * eta[x] * nu[x] * t ^ (nu[x] - 1) * phi_x(x, k, L)
  }

  # Summed cumulative hazard
  sum_cum_haz <- function(u, t, k, L) {
    sum(sapply(x, function(x) {
      at_risk(x - 1, k, L) * eta[x] * phi_x(x, k, L) * ((t + u) ^ nu[x] - t ^ nu[x])
    }))}


  # Inverse summed cumulative hazard function
  inverse_sc_haz <- function(p, t, k, L, lower_bound = 10^-15, upper_bound = 100) {
    root_function <- function(u) sum_cum_haz(u, t, k, L) - p
    stats::uniroot(root_function, lower = lower_bound, upper = upper_bound)$root
  }

  # Event probabilities
  probs <- function(t, k, L){
    probs <- sapply(x, function(x) lambda_x(x, t, k, L))
    summ <- sum(probs)
    probs / summ
  }

  # Simulation
  res <- data.table(ID = numeric(), Time = numeric(), Delta = numeric(), L0 = numeric(),
                    L = numeric(), A = numeric())

  for(j in 1:N){
    # Draw
    L0 <- stats::runif(1)
    A <- stats::rbinom(1, 1, 0.5)

    # Initialize
    T_k1 <- 0
    Delta <- -1
    L <- numeric()
    L[1] <- 0

    # Iterate
    Ts <- c()
    Deltas <- c()
    # Indicator for number of events
    k <- 0

    while(! Delta %in% term_deltas){
      # Simulate time
      V <- stats::runif(1)
      W_k1 <- inverse_sc_haz(-log(V), T_k1, k, L[k + 1])
      T_k1 <- T_k1 + W_k1
      Ts[k + 1] <- T_k1
      # Simulate event
      Deltas[k + 1] <- sample(x, size = 1, prob = probs(T_k1, k, L[k + 1])) - 1
      Delta <- Deltas[k + 1]
      # Update number of events and covariate change indicator
      k <- k + 1
      if(Delta == 3) L[k + 1] <- 1
      else L[k + 1] <- L[k]
    }
    jth_res <- data.table(ID = rep(j, k),
                          Time = Ts,
                          Delta = Deltas,
                          L0 = rep(L0,k),
                          L = L[1:k],
                          A = rep(A,k))
    res <- rbind(res, jth_res)
  }
  return(res)
}
