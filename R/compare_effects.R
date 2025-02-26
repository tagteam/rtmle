#' `compare_effects` is a function for computing a user specified estimate for three
#' different effects of A0 on the intensity of L, and for a varying parameter.
#'
#' @title Compute estimate for varying effects
#'
#' @param estimator A user specified estimator. A function that takes as input data and N,
#' and calculates an estimate.
#' @param N Number of individuals simulated
#' @param beta_L0_L Parameter for the effect of L0 on L.
#' @param beta_L_D Parameter for the effect of L on D.
#' @param beta_A0_L Parameter for the effect of L0 on L.
#' @param beta_L0_D Parameter for the effect of L0 on D.
#' @param nu Vector of length 4 of scale parameters for the Weibull intensity with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 1.1 for all events.
#' @param eta Vector of length 4 of shape parameters for the Weibull intensity. Default is set to
#' 0.1 for all events.
#' @param cens Indicator for whether censoring is present
#' @param setting Setting specifier
#'
#' @return Plot
#' @export
#'
#' @examples
#' # First we define an estimator
#' estimator1 <- function(data, N) {
#'
#' T2D_events <- data[Delta == 3]
#' return(c(nrow(T2D_events[A0 == 1 & Time < 5])/N, nrow(T2D_events[A0 == 0 & Time < 5])/N))
#'
#' }
#'
#' res1 <- compare_effects(estimator = estimator1, N = 100, beta_L_D = seq(0,1,by = 0.1))

compare_effects <- function(estimator,
                            N = 1000,
                            beta_L0_L = 2,
                            beta_L_D = 1,
                            beta_A0_L = 0,
                            beta_L0_D = 0,
                            nu = rep(1.1,4),
                            cens = 0,
                            eta = rep(0.1,4),
                            setting = 1) {

  B <-  max(length(beta_L0_L), length(beta_L_D), length(beta_A0_L), length(beta_L0_D))

  if(length(beta_L0_L) < B) beta_L0_L <- rep(beta_L0_L, B)
  if(length(beta_L_D) < B) beta_L_D <- rep(beta_L_D, B)
  if(length(beta_A0_L) < B) beta_A0_L <- rep(beta_A0_L, B)
  if(length(beta_L0_D) < B) beta_L0_D <- rep(beta_L0_D, B)

  if(setting == 1){
    for(b in 1:B){
      data_a <- sim_data_setting2(N, beta_L_D = beta_L_D[b], beta_A0_L = beta_A0_L[b],
                                 beta_L0_L = beta_L0_L[b], beta_A0_D = 0, nu = nu,
                                 cens = 0, eta = eta, beta_L0_D = beta_L0_D[b])
      data_b <- sim_data_setting2(N, beta_L_D = beta_L_D[b], beta_A0_L = beta_A0_L[b],
                                   beta_L0_L = beta_L0_L[b], beta_A0_D = -0.1, nu = nu,
                                   cens = 0, eta = eta, beta_L0_D = beta_L0_D[b])
      data_c <- sim_data_setting2(N, beta_L_D = beta_L_D[b], beta_A0_L = beta_A0_L[b],
                                   beta_L0_L = beta_L0_L[b], beta_A0_D = -0.2, nu = nu,
                                   cens = 0, eta = eta, beta_L0_D = beta_L0_D[b])

      if(b == 1){
        est0 <- estimator(data_a, N)
        est_length <- length(est0)
        estimates <- matrix(nrow = B, ncol = 3 * est_length)
      }

      estimates[b, 1 : est_length] <- estimator(data_a, N)
      estimates[b, (est_length + 1) : (2 * est_length)] <- estimator(data_b, N)
      estimates[b, (2 * est_length + 1) : (3 * est_length)] <- estimator(data_c, N)
    }
  }
  else{
    for(b in 1:B){
      data_a <- sim_data_setting2(N, beta_L_D = 0, beta_A0_L = beta_A0_L[b],
                                 beta_L0_L = beta_L0_L[b], beta_A0_D = 0, nu = nu,
                                 cens = 0, eta = eta, beta_L0_D = beta_L0_D[b])
      data_b <- sim_data_setting2(N, beta_L_D = 0.5, beta_A0_L = beta_A0_L[b],
                                   beta_L0_L = beta_L0_L[b], beta_A0_D = 0, nu = nu,
                                   cens = 0, eta = eta, beta_L0_D = beta_L0_D[b])
      data_c <- sim_data_setting2(N, beta_L_D = 1, beta_A0_L = beta_A0_L[b],
                                   beta_L0_L = beta_L0_L[b], beta_A0_D = 0, nu = nu,
                                   cens = 0, eta = eta, beta_L0_D = beta_L0_D[b])

      if(b == 1){
        est0 <- estimator(data_a, N)
        est_length <- length(est0)
        estimates <- matrix(nrow = B, ncol = 3 * est_length)
      }

      estimates[b, 1 : est_length] <- estimator(data_a, N)
      estimates[b, (est_length + 1) : (2 * est_length)] <- estimator(data_b, N)
      estimates[b, (2 * est_length + 1) : (3 * est_length)] <- estimator(data_c, N)
    }
  }

  return(estimates)
}
