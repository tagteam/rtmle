### plot.rtmle.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (16:42) 
## Version: 
## Last-Updated: maj 28 2026 (12:26) 
##           By: Thomas Alexander Gerds
##     Update #: 68
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
#' @param targets Character vector selecting target labels from the
#'   \code{Target} column. If \code{NULL}, all stored targets are plotted.
#' @param xlim Limits for the x-axis.
#' @param ylim Limits for the y-axis.
#' @param y_breaks Breaks for the y-axis.
#' @param x_breaks Breaks for the x-axis.
#' @param position_atrisk Vector of positions on the x-axis where numbers at-risk are shown below the graph.
#' @param conf_int Logical. If \code{TRUE}, add confidence bands to the graph.
#' @param bootstrap_conf_int Logical. If \code{TRUE}, add shaded confidence
#'   bands from the cheap-bootstrap results stored in
#'   \code{object$estimate$Cheap_bootstrap}. The selected analysis must have
#'   corresponding bootstrap results.
#' @param ... Not used.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{run_rtmle}}, \code{\link{summary.rtmle}},
#'   \code{\link{plot_adherence}}, \code{\link{plot_IPW}}
#' @examples
#' data(rtmle_object)
#' p <- ggplot2::autoplot(rtmle_object, xlim = c(0, 4))
#' class(p)
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs
#'     theme_minimal facet_wrap
#' @rdname plot.rtmle
#' @method autoplot rtmle
#' @export
autoplot.rtmle <- function(object,
                           analysis = NULL,
                           targets = NULL,
                           xlim,
                           ylim,
                           y_breaks,
                           x_breaks,
                           position_atrisk,
                           conf_int,
                           bootstrap_conf_int = FALSE,
                           ...) {
    B=Bootstrap_lower=Bootstrap_upper=.plot_group <- NULL
    Estimate=Lower=N=Regime=Target=Time_horizon=Upper=last_interval <- NULL
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#F0E442")
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' must be installed to use autoplot.rtmle().")
    }
    time_unit <- if (is.null(object$time_unit)) "Time" else object$time_unit

    if (length(bootstrap_conf_int) != 1L ||
        !is.logical(bootstrap_conf_int) ||
        is.na(bootstrap_conf_int)) {
        stop("bootstrap_conf_int must be a single TRUE/FALSE value.")
    }

    # Collect the requested analysis into a data.table.
    analysis_name <- if (is.null(analysis)) "Main_analysis" else analysis
    if (length(analysis_name) != 1L ||
        is.na(analysis_name) ||
        !analysis_name %in% names(object$estimate)) {
        stop("Analysis '", paste(analysis_name, collapse = ", "),
             "' not found in object$estimate")
    }
    est <- object$estimate[[analysis_name]]
    if (!(inherits(est, "data.frame") || is.matrix(est))) {
        stop("Analysis '", analysis_name,
             "' does not contain a tabular set of estimates.")
    }
    est <- data.table::as.data.table(data.table::copy(est))
    required_estimate_columns <- c(
        "Regime", "Time_horizon", "Estimate", "Lower", "Upper"
    )
    missing_estimate_columns <- setdiff(required_estimate_columns, names(est))
    if (length(missing_estimate_columns) > 0L) {
        stop("Analysis '", analysis_name, "' is missing required column(s): ",
             paste(missing_estimate_columns, collapse = ", "))
    }
    if (!("Target" %in% names(est))) {
        est[, Target := ""]
    }
    if (!is.null(targets)) {
        if (length(targets) == 0L) {
            stop("targets must contain at least one target label or be NULL.")
        }
        target_labels <- as.character(targets)
        available_targets <- unique(as.character(est$Target))
        missing_targets <- setdiff(target_labels, available_targets)
        if (length(missing_targets) > 0L) {
            stop("Target(s) not found in analysis '", analysis_name, "': ",
                 paste(missing_targets, collapse = ", "),
                 ". Available targets: ",
                 paste(available_targets, collapse = ", "), ".")
        }
        est <- est[as.character(Target) %in% target_labels]
    }
    if (NROW(est) == 0L) {
        stop("No estimates remain after filtering by targets.")
    }

    bootstrap_est <- NULL
    if (isTRUE(bootstrap_conf_int)) {
        bootstrap_results <- object$estimate[["Cheap_bootstrap"]]
        if (is.null(bootstrap_results) ||
            !is.list(bootstrap_results) ||
            is.null(bootstrap_results[[analysis_name]])) {
            stop("Bootstrap confidence intervals were requested, but no ",
                 "Cheap_bootstrap results were found for analysis '",
                 analysis_name, "'.")
        }
        bootstrap_est <- data.table::as.data.table(data.table::copy(
            bootstrap_results[[analysis_name]]
        ))
        required_bootstrap_columns <- c(
            "Regime", "Time_horizon", "Bootstrap_lower", "Bootstrap_upper"
        )
        missing_bootstrap_columns <- setdiff(
            required_bootstrap_columns,
            names(bootstrap_est)
        )
        if (length(missing_bootstrap_columns) > 0L) {
            stop("Cheap_bootstrap results for analysis '", analysis_name,
                 "' are missing required column(s): ",
                 paste(missing_bootstrap_columns, collapse = ", "), ".")
        }
        if (!("Target" %in% names(bootstrap_est))) {
            bootstrap_est[, Target := ""]
        }
        if (!is.null(targets)) {
            bootstrap_est <- bootstrap_est[
                as.character(Target) %in% as.character(targets)
            ]
        }
        if (NROW(bootstrap_est) == 0L) {
            stop("No bootstrap confidence intervals remain after filtering ",
                 "by targets.")
        }
        if ("B" %in% names(bootstrap_est)) {
            bootstrap_est <- bootstrap_est[
                , .SD[which.max(B)],
                by = c("Target", "Regime", "Time_horizon")
            ]
        }
        bootstrap_est[, .plot_group := interaction(
            as.character(Target), as.character(Regime), drop = TRUE
        )]
        bootstrap_zero <- unique(bootstrap_est[, list(Target, Regime)])
        bootstrap_zero[, c(
            "Time_horizon", "Bootstrap_lower", "Bootstrap_upper"
        ) := list(0, 0, 0)]
        bootstrap_zero[, .plot_group := interaction(
            as.character(Target), as.character(Regime), drop = TRUE
        )]
        bootstrap_est <- data.table::rbindlist(
            list(
                bootstrap_zero,
                bootstrap_est[, list(
                    Target,
                    Regime,
                    Time_horizon,
                    Bootstrap_lower,
                    Bootstrap_upper,
                    .plot_group
                )]
            ),
            use.names = TRUE,
            fill = TRUE
        )
    }

    # Assume that estimates are zero at time zero.
    est[, .plot_group := interaction(
        as.character(Target), as.character(Regime), drop = TRUE
    )]
    est_zero <- unique(est[, list(Target, Regime)])
    est_zero[, c("Time_horizon", "Estimate", "Lower", "Upper") := list(0, 0, 0, 0)]
    est_zero[, .plot_group := interaction(
        as.character(Target), as.character(Regime), drop = TRUE
    )]
    est <- data.table::rbindlist(
        list(
            est_zero,
            est[, list(
                Target,
                Regime,
                Time_horizon,
                Estimate,
                Lower,
                Upper,
                .plot_group
            )]
        ),
        use.names = TRUE,
        fill = TRUE
    )
    # Build ggplot
    p <- ggplot2::ggplot(est, ggplot2::aes(x = Time_horizon,
                                           y = Estimate,
                                           color = Regime,
                                           fill = Regime,
                                           group = .plot_group))
    if (missing(x_breaks)) {
        x_breaks <- object$time_grid
    }
    # getting data for numbers at-risk below the graph
    if (missing(position_atrisk)){
        position_atrisk <- x_breaks
    }
    atrisk <- object$followup[,.N,keyby = last_interval]
    atrisk[,N := NROW(object$followup)-cumsum(c(0,N[-length(N)]))]
    atrisk <- atrisk[last_interval%in%position_atrisk]
    if (match("Lower",names(est),nomatch = 0) >0 &&
        ((missing(conf_int) ||
          (length(conf_int)>0 && conf_int != FALSE)))){
        p <- p + ggplot2::geom_ribbon(
            ggplot2::aes(ymin = Lower, ymax = Upper),
            linetype = 0,
            alpha = 0.2
        )
    }
    if (isTRUE(bootstrap_conf_int)) {
        p <- p + ggplot2::geom_ribbon(
            data = bootstrap_est,
            mapping = ggplot2::aes(
                x = Time_horizon,
                ymin = Bootstrap_lower,
                ymax = Bootstrap_upper,
                fill = Regime,
                group = .plot_group
            ),
            inherit.aes = FALSE,
            linetype = 0,
            alpha = 0.12,
            show.legend = FALSE
        )
    }
    p <- p+ggplot2::geom_line()+ggplot2::geom_point()
    p <- p+ ggplot2::labs(x = time_unit,
                          y = "Estimated risk",
                          title = "")
    p <- p+ggplot2::scale_fill_manual(values = cbbPalette)
    p <- p+ggplot2::scale_color_manual(values = cbbPalette)
    p <- p+ggplot2::theme_minimal(base_size = 12)
    # axes
    if (missing(ylim)) ylim <- c(0,1)
    if (missing(xlim)) xlim <- range(x_breaks, na.rm = TRUE)
    ## p <- p+ggplot2::ylim(0,1)
    if (missing(y_breaks)) y_breaks <- seq(ylim[1],ylim[2],abs(ylim[2]-ylim[1])/4)
    p <- p+ggplot2::scale_y_continuous(limits = ylim,breaks = y_breaks,labels = paste0(100*y_breaks,"%"))
    p <- p+ggplot2::scale_x_continuous(breaks = x_breaks, labels = object$time_grid_labels[x_breaks+1])
    p <- p+ggplot2::coord_cartesian(ylim = ylim,xlim = xlim,clip = 'off')
    # FIXME: this should not be necessary but at some point we want to add the
    #        number of people who actually follow the regime and are at-risk
    atrisk[,Regime := ""]
    p <- p+ggplot2::geom_text(data = atrisk,
                              mapping = ggplot2::aes(x = last_interval,
                                                     y = I(0),
                                                     vjust = 9,
                                                     label = N,
                                                     fill = NULL,
                                                     colour = NULL),
                              inherit.aes = FALSE,
                              show.legend = FALSE)
    # space for atrisk data
    p <- p+ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,6,1), "lines"))
    p <- p+ggplot2::ylab("Risk")+ggplot2::xlab(time_unit)
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
plot.rtmle <- function(x,
                       analysis = NULL,
                       targets = NULL,
                       bootstrap_conf_int = FALSE,
                       ...) {
    autoplot.rtmle(
        x,
        analysis = analysis,
        targets = targets,
        bootstrap_conf_int = bootstrap_conf_int,
        ...
    )
}


######################################################################
### plot.rtmle.R ends here
