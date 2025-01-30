#' `plot_compare` is a function for plotting an estimates. It takes a matrix of estimates, and plots.
#'
#' @title Plot comparison of estimates.
#'
#' @param estimates A matrix B times 6 matrix with estimates data where A0 = 0 or 1
#' and where the effect of the drug varies from no effect, small effect to large effect.
#' @param plot_no For now the possibility is 1 or 2. Plot_no = 1 lets us plot the different
#' estimates. While plot_no = 2, lets us plot differences.
#' @param diff_betas The values of the varying beta coefficient.
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

plot_compare <- function(estimates, plot_no = 1, diff_betas = seq(0,1,by = 0.1)) {

  ests <- group <- No_effect <- Effect <- Large_effect <- NULL

  if(plot_no == 1){
    plot_data <- data.table(ests = c(estimates),
                            group = c(rep("A0 = 0, No effect", 11), rep("A0 = 1, No effect", 11),
                                      rep("A0 = 0, Small effect", 11), rep("A0 = 1, Small effect", 11),
                                      rep("A0 = 0, Large effect", 11), rep("A0 = 1, Large effect", 11)),
                            beta = rep(diff_betas,6))

    p <-  ggplot2::ggplot(plot_data)+
      ggplot2::geom_line(ggplot2::aes(x = beta, y = ests, colour = group))
  }

  else if(plot_no == 2){
    plot_data <- data.table(No_effect = estimates[,2] - estimates[,1],
               Effect = estimates[,4] - estimates[,3],
               Large_effect = estimates[,6] - estimates[,5],
               beta = rep(diff_betas,2))

      p <-  ggplot2::ggplot(plot_data)+
        ggplot2::geom_line(ggplot2::aes(x = beta, y = No_effect, colour = "No effect of drug on T2D"))+
        ggplot2::geom_line(ggplot2::aes(x = beta, y = Effect, colour = "Small effect of drug on T2D"))+
        ggplot2::geom_line(ggplot2::aes(x = beta, y = Large_effect, colour = "Large effect of drug on T2D"))+
        ggplot2::scale_color_manual(values = c("darkolivegreen", "darkblue", "hotpink"))
  }
  return(p)
}
