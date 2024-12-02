#' Vectorized function to simulate event data. The number of events simulated corresponds
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
#' #vec_sim_event_data2(N = 10)
#'
#' #Not working
#' #Work in progress

vec_sim_event_data2 <- function(N,                  # Number of individuals
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
    at_risk <- function(x) {
      # If you have not died yet or been censored yet, you are at risk for dying or being censored
      if(x == 1 | x == 2) return(rep(1, N))
      # You are only at risk for an operation if you have not had an operation yet
      else if(x == 0) return(as.numeric(k == 0 | (k == 1 & L == 1)))
      # You are only at risk for a change in the covariate process if you have not had a change yet
      else return(as.numeric(L == 0))
    }
  }


  # Events
  x <- 1:ncol(beta)

  # Intensities
  phi_x <- function(x) {
    exp(L0 * beta[1,x] + k * beta[2,x] + A * beta[3,x] + L * beta[4,x])
  }

  lambda_x <- function(x, t) {
    at_risk(x - 1) * eta[x] * nu[x] * t ^ (nu[x] - 1) * phi_x(x)
  }

  # Summed cumulative hazard
  sum_cum_haz <- function(u, t) {
    sum(sapply(x, function(x) {
      at_risk(x - 1) * eta[x] * phi_x(x) * ((t + u) ^ nu[x] - t ^ nu[x])
    }))}


  # Inverse summed cumulative hazard function
  inverse_sc_haz <- function(p, t, lower_bound = 10^-15, upper_bound = 100) {
    root_function <- function(u) sum_cum_haz(u, t) - p
    stats::uniroot(root_function, lower = lower_bound, upper = upper_bound)$root
  }

  # Event probabilities
  probs <- function(t){
    probs <- sapply(x, function(x) lambda_x(x, t))
    summ <- sum(probs)
    probs / summ
  }

  # Draw
  L0 <- stats::runif(N)
  A <- stats::rbinom(N, 1, 0.5)

  # Initialize
  T_k <- rep(0,N)
  Delta <- -1
  L <- rep(0, N)
  k <- rep(0,N)
  alive <- 1:N

  res <- data.table()

  while(length(alive) != 0){
    # Simulate time
    V <- stats::runif(N)
    W <- inverse_sc_haz(-log(V), T_k)
    T_k[alive] <- T_k[alive] + W

    # Simulate event
    Deltas <- sapply(alive, function(i) sample(x, size = 1, prob = probs(T_k[i], i)) - 1)

    kth_event <- data.table(ID = alive,
                            Time = T_k[alive],
                            Delta = Deltas,
                            L0 = L0[alive],
                            L = L[alive],
                            A = A[alive])

    res <- rbind(res, kth_event)

    # Update number of events and covariate change
    k <- k + 1
    L[alive][Deltas == 3] <- 1

    # Who is still alive and uncensored?
    alive <- alive[! Deltas %in% term_deltas]
  }

  ID <- NULL

  return(res[order(ID)])
}
