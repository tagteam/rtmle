#' `plot_compare` is a function for plotting an estimates. It takes a matrix of estimates, and plots.
#'
#' @title Plot comparison of estimates.
#'
#' @param estimates A matrix with estimates. Columns 1,3,... are for A0 = 0 and
#' columns 2,4,... for  A0 = 1
#' @param diff_betas The values of the varying beta coefficient.
#' @param lower A matrix with lower confidence bands. Columns 1,3,... are for A0 = 0 and
#' columns 2,4,... for  A0 = 1
#' @param upper A matrix with lower confidence bands. Columns 1,3,... are for A0 = 0 and
#' columns 2,4,... for  A0 = 1
#' @param groups The strata for the different estimates
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
#' plot_compare(res1, diff_betas = seq(0,1,by = 0.1))

plot_compare <- function(estimates, diff_betas = seq(0,1,by = 0.1),
                         lower = NULL,
                         upper = NULL,
                         groups = c(0, -0.1, -0.2)) {

  ests <- A0 <- NULL

  B <- length(diff_betas)
  n_cols <- ncol(estimates)

  plot_data <- data.table(ests = c(estimates),
                          A0 = rep(c(rep("0", B), rep("1", B)), n_cols / 2),
                          beta = rep(diff_betas, n_cols),
                          group = c(rep(groups[1], 2 * B), rep(groups[2], 2* B), rep(groups[3], 2* B)))

  if (!is.null(upper) && !is.null(lower)) {
    plot_data[, upper := c(upper)]
    plot_data[, lower := c(lower)]
  }

  p <-  ggplot2::ggplot(plot_data)+
    ggplot2::geom_line(ggplot2::aes(x = beta, y = ests, colour = A0))+
    ggplot2::facet_wrap(~ group)

  # Add shading only if both upper and lower are provided
  if (!is.null(upper) && !is.null(lower)) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(x = beta, ymin = lower, ymax = upper,
                                        fill = A0), alpha = 0.2)
  }

  return(p)
}
