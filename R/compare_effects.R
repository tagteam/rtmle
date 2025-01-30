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
#' @param beta_A0_D Parameter for the effect of L0 on L. One of the three parameters is
#' allowed to vary.
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

compare_effects <- function(estimator, N = 1000, beta_L0_L = 1, beta_L_D = 1, beta_A0_D = 0) {

  B <-  max(length(beta_L0_L), length(beta_L_D), length(beta_A0_D))

  if(length(beta_L0_L) < B) beta_L0_L <- rep(beta_L0_L, B)
  if(length(beta_L_D) < B) beta_L_D <- rep(beta_L_D, B)
  if(length(beta_A0_D) < B) beta_A0_D <- rep(beta_A0_D, B)

  estimates <- matrix(nrow = B, ncol = 2*3)

  for(b in 1:B){
    print(b)
    data0 <- sim_data_setting2(N, beta_L_D = beta_L_D[b], beta_A0_D = beta_A0_D[b],
                               beta_L0_L = beta_L0_L[b], beta_A0_L = 0)
    data0.5 <- sim_data_setting2(N, beta_L_D = beta_L_D[b], beta_A0_D = beta_A0_D[b],
                                 beta_L0_L = beta_L0_L[b], beta_A0_L = -0.5)
    data1 <- sim_data_setting2(N, beta_L_D = beta_L_D[b], beta_A0_D = beta_A0_D[b],
                               beta_L0_L = beta_L0_L[b], beta_A0_L = -1)

    estimates[b,c(1,2)] <- estimator(data0, N)
    estimates[b,c(3,4)] <- estimator(data0.5, N)
    estimates[b,c(5,6)] <- estimator(data1, N)
  }
  return(estimates)
}
