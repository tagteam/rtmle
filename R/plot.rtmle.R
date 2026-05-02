### plot.rtmle.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (16:42) 
## Version: 
## Last-Updated: apr 25 2026 (06:56) 
##           By: Thomas Alexander Gerds
##     Update #: 35
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Plot rtmle objects
#' @description Plot risk estimates and risk contrasts from fitted
#'   \code{rtmle} objects.
#' @name plot.rtmle
#' @aliases autoplot.rtmle
#' @param object A prepared and fitted object of class \code{"rtmle"}.
#' @param x A prepared and fitted object of class \code{"rtmle"}.
#' @param analysis Name of the analysis. If \code{NULL}, use
#'   \code{object$estimate[["Main_analysis"]]}.
#' @param xlim Limits for the x-axis.
#' @param ylim Limits for the y-axis.
#' @param y_breaks Breaks for the y-axis.
#' @param x_breaks Breaks for the x-axis.
#' @param position_atrisk Vector of positions on the x-axis where numbers at-risk are shown below the graph.
#' @param conf_int Logical. If \code{TRUE}, add confidence bands to the graph.
#' @param ... Not used.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{run_rtmle}}, \code{\link{summary.rtmle}},
#'   \code{\link{plot_adherence}}, \code{\link{plot_IPW}}
#' @examples
#' x <- list(
#'   estimate = list(Main_analysis = data.table::data.table(
#'     Time_horizon = c(1, 2, 1, 2),
#'     Time = c(1, 2, 1, 2),
#'     Protocol = rep(c("Always_A", "Never_A"), each = 2),
#'     Estimate = c(.12, .20, .18, .30),
#'     Lower = c(.08, .15, .12, .24),
#'     Upper = c(.18, .27, .25, .38))),
#'   followup = data.table::data.table(id = 1:5,
#'                                     last_interval = c(0, 1, 2, 2, 2)),
#'   time_grid_scale = 0:2)
#' class(x) <- "rtmle"
#' p <- ggplot2::autoplot(x, xlim = c(0, 2))
#' class(p)
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs
#'     theme_minimal facet_wrap
#' @rdname plot.rtmle
#' @method autoplot rtmle
#' @export
autoplot.rtmle <- function(object,
                           analysis = NULL,
                           xlim,
                           ylim,
                           y_breaks,
                           x_breaks,
                           position_atrisk,
                           conf_int,
                           ...) {
    Estimate=Lower=N=Protocol=Time_horizon=Upper=last_interval <- NULL
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#F0E442")
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' must be installed to use autoplot.rtmle().")
    }
    ## if (!requireNamespace("pammtools", quietly = TRUE)) {
        ## stop("Package 'pammtools' must be installed to use autoplot.rtmle().")
    ## }

    # Collect all analyses into one data.frame
    if (!is.null(analysis)) {
        if (!analysis %in% names(object$estimate)) {
            stop("Analysis '", analysis, "' not found in object$estimate")
        }
        est <- object$estimate[[analysis]]
    }else{
        est <- object$estimate[["Main_analysis"]] 
    }
    # Build ggplot
    p <- ggplot2::ggplot(est, ggplot2::aes(x = Time_horizon,
                                           y = Estimate,
                                           color = Protocol,
                                           fill = Protocol,
                                           group = Protocol))
    # getting data for numbers at-risk below the graph
    if (missing(position_atrisk)){
        position_atrisk <- object$time_grid_scale
    }
    atrisk_times <- data.table::data.table(time = position_atrisk)
    atrisk <- object$followup[,.N,keyby = last_interval]
    atrisk[,N := NROW(object$followup)-cumsum(c(0,N[-length(N)]))]
    ## p <- p+ggplot2::geom_step()
    p <- p+ggplot2::geom_line()+ggplot2::geom_point()
    color_variable <- "Protocol"
    if (match("Lower",names(est),nomatch = 0) >0 &&
        ((missing(conf_int) ||
          (length(conf_int)>0 && conf_int != FALSE)))){
        ## p <- p+ pammtools::geom_stepribbon(ggplot2::aes(ymin = Lower,ymax = Upper,fill = !! if (length(color_variable)>0) {rlang::sym(color_variable)}else{NULL},colour = !! if (length(color_variable)>0) {rlang::sym(color_variable)}else{NULL}),linetype = 0,alpha = 0.2)
        p <- p+ geom_ribbon(ggplot2::aes(ymin = Lower,ymax = Upper,fill = !! if (length(color_variable)>0) {rlang::sym(color_variable)}else{NULL},colour = !! if (length(color_variable)>0) {rlang::sym(color_variable)}else{NULL}),linetype = 0,alpha = 0.2)
    }
    p <- p+ ggplot2::labs(x = "Time horizon",
                          y = "Estimated risk",
                          title = "")
    p <- p+ggplot2::scale_fill_manual(values = cbbPalette)
    p <- p+ggplot2::scale_color_manual(values = cbbPalette)
    p <- p+ggplot2::theme_minimal(base_size = 12)
    # axes
    if (missing(ylim)) ylim <- c(0,1)
    if (missing(xlim)) xlim <- c(0,max(est$Time))
    ## p <- p+ggplot2::ylim(0,1)
    if (missing(y_breaks)) y_breaks <- seq(ylim[1],ylim[2],abs(ylim[2]-ylim[1])/4)
    p <- p+ggplot2::scale_y_continuous(limits = ylim,breaks = y_breaks,labels = paste0(100*y_breaks,"%"))
    p <- p+ggplot2::coord_cartesian(ylim = ylim,xlim = xlim,clip = 'off')
    # FIXME: this should not be necessary but at some point we want to add the
    #        number of people who actually follow the regime and are at-risk
    atrisk[,Protocol := ""]
    p <- p+ggplot2::geom_text(data = atrisk,
                              mapping = ggplot2::aes(x = last_interval,
                                                     y = I(0),
                                                     vjust = 9,
                                                     label = N,
                                                     fill = NULL,
                                                     colour = NULL),
                              show.legend = FALSE)
    # space for atrisk data
    p <- p+ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,6,1), "lines"))
    p <- p+ggplot2::ylab("Risk")+ggplot2::xlab("Time")
    p <- p+ggplot2::annotate("text",
                             x = 0,
                             y = I(0),
                             vjust = ggplot2::unit(7, "lines"),
                             label = "Number at risk")
    return(p)
}
#' @method plot rtmle
#' @rdname plot.rtmle
#' @export
plot.rtmle <- function(x, ...) {
  print(autoplot.rtmle(x, ...))
}


######################################################################
### plot.rtmle.R ends here
