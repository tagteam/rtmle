#' Plot estimates over time with one row per stratum
#' and one column per target.
#'
#' @param object A fitted \code{stratified_rtmle} object.
#' @param targets Character vector of target names. The default uses all fitted
#'   targets.
#' @param target_parameter The quantity to plot. Choose \code{"Risk"},
#'   \code{"Risk_difference"}, or \code{"Risk_ratio"}.
#' @param reference Optional named list of reference regimes, one per target,
#'   as in \code{\link{summary.rtmle}}.
#' @param conf_int Logical. Draw pointwise 95\% confidence bands if \code{TRUE}.
#' @param xlim,ylim Optional limits for the horizontal and vertical axes.
#'   Horizontal limits are specified as time-grid indices. Vertical limits are
#'   specified on the scale of the estimate, for example \code{c(0, 0.2)} for
#'   risks between 0 and 20\%.
#' @param ... Additional arguments such as \code{digits} passed to
#'   \code{\link{summary.rtmle}} for each fit. In the plot method, arguments
#'   are passed to \code{autoplot()}.
#' @return returns a
#'   \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{stratify.rtmle}}, \code{\link{summary.rtmle}}
#' @examples
#' data(rtmle_object)
#' by_age <- stratify.rtmle(rtmle_object, age, breaks = c(-Inf, 50, Inf),
#'                         time_horizon = 1:2)
#' plot(by_age)
#'
#' # Compare Always_A with Always_B within each age group for one target.
#' target_name <- "Treatment contrast"
#' reference <- list("Treatment contrast" = "Always_B")
#' plot(by_age, targets = target_name, reference = reference,
#'      target_parameter = "Risk_difference")
#' @name autoplot.stratified_rtmle
NULL

#' @rdname autoplot.stratified_rtmle
#' @export
autoplot.stratified_rtmle <- function(object,
                                     targets = NULL,
                                     target_parameter = c("Risk", "Risk_difference", "Risk_ratio"),
                                     reference = NULL,
                                     conf_int = TRUE,
                                     xlim = NULL,
                                     ylim = NULL,
                                     ...) {
    Estimate <- Lower <- Regime <- Series <- Stratum <- Target <- Time_horizon <- Upper <- NULL
    target_parameter <- match.arg(target_parameter)
    estimates <- if (is.null(targets)) {
        summary(object, reference = reference, ...)
    } else {
        summary(object, targets = targets, reference = reference, ...)
    }
    data.table::setnames(estimates, 1L, "Stratum")
    estimates <- estimates[estimates$Target_parameter == target_parameter]
    if (nrow(estimates) == 0L) {
        stop("No estimates are available for '", target_parameter, "'.")
    }
    if (target_parameter == "Risk") {
        # As in plot.rtmle(), the risk starts at zero at baseline.
        baseline <- unique(estimates[, list(Stratum, Target, Regime)])
        data.table::set(baseline, j = c("Time_horizon", "Estimate", "Lower", "Upper"),
                        value = list(0, 0, 0, 0))
        estimates <- data.table::rbindlist(list(baseline, estimates),
                                          use.names = TRUE, fill = TRUE)
        data.table::set(estimates, j = "Series", value = estimates$Regime)
        legend_title <- "Regime"
    } else {
        data.table::set(estimates, j = "Series",
                        value = paste(estimates$Regime, "vs", estimates$Reference))
        legend_title <- "Regime contrast"
    }
    fit <- object[[1L]]
    time_breaks <- sort(unique(estimates$Time_horizon))
    time_labels <- if (is.null(fit$time_grid_labels)) {
        time_breaks
    } else {
        fit$time_grid_labels[match(time_breaks, fit$time_grid)]
    }
    time_unit <- if (is.null(fit$time_unit)) "Time" else fit$time_unit
    variable <- attr(object, "variable", exact = TRUE)
    stratum_labels <- function(labels) {
        vapply(paste(variable, "=", labels), function(label) {
            paste(strwrap(label, width = 28L), collapse = "\n")
        }, character(1L), USE.NAMES = FALSE)
    }
    p <- ggplot2::ggplot(estimates, ggplot2::aes(Time_horizon, Estimate,
                                                 colour = Series, fill = Series))
    if (target_parameter != "Risk") {
        p <- p + ggplot2::geom_hline(yintercept = as.numeric(target_parameter == "Risk_ratio"),
                                    linetype = 2L, colour = "grey50")
    }
    if (isTRUE(conf_int)) {
        limits <- ggplot2::aes(ymin = Lower, ymax = Upper)
        interval <- if (length(time_breaks) > 1L) {
            ggplot2::geom_ribbon(limits, alpha = 0.2, colour = NA, na.rm = TRUE)
        } else {
            ggplot2::geom_linerange(limits, na.rm = TRUE)
        }
        p <- p + interval
    }
    if (length(time_breaks) > 1L) {
        p <- p + ggplot2::geom_line(na.rm = TRUE)
    }
    p <- p + ggplot2::geom_point(na.rm = TRUE) +
        ggplot2::facet_grid(rows = ggplot2::vars(Stratum),
                            cols = ggplot2::vars(Target),
                            labeller = ggplot2::labeller(
                                Stratum = stratum_labels,
                                .default = ggplot2::label_wrap_gen(width = 28L))) +
        ggplot2::scale_x_continuous(breaks = time_breaks, labels = time_labels) +
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
        ggplot2::labs(x = time_unit, colour = legend_title, fill = legend_title) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
    switch(target_parameter,
           Risk = p + ggplot2::scale_y_continuous("Risk", labels = function(x) paste0(100 * x, "%")),
           Risk_difference = p + ggplot2::scale_y_continuous("Risk difference (percentage points)",
                                                            labels = function(x) 100 * x),
           Risk_ratio = p + ggplot2::ylab("Risk ratio"))
}

