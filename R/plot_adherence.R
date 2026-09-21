### plot_adherence.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: dec 11 2025 (10:23) 
## Version: 
## Last-Updated: sep 18 2026 (11:39)
##           By: Thomas Alexander Gerds
##     Update #: 41
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Plot cumulative non-adherence by regime
#'
#' Uses \code{\link{adherence}} to collect time to first deviation from a
#' treatment regime for each regime in a fitted \code{rtmle} object
#' \code{x}, allowing for right censoring and competing risks (e.g.,
#' death/outcome). The result is plotted as the cumulative incidence (in
#' percent) of non-adherence over follow-up time, stratified by regime.
#'
#' @param x An object containing regime-specific adherence information and follow-up
#'   data. Must include at least \code{x$regimes} (a named list where each element has
#'   \code{$intervention_match}), \code{x$followup} (with
#'   \code{last_interval}), \code{x$prepared_data} (optional; used for
#'   censoring indicators), and \code{x$names$censoring}.
#' @param regimes Optional names of the regimes to plot. If omitted, use
#'   all regimes that already have an \code{intervention_match} table.
#' @param ... Currently unused. Included for future extensions.
#'
#' @details
#' The underlying \code{\link{adherence}} function restricts each regime to initiators (those with
#' \code{intervention_match[,1] == 1}). It then identifies:
#' \itemize{
#'   \item \code{first_deviation}: the first interval where \code{intervention_match} equals 0.
#'   \item \code{censored_time}: the first interval marked \code{"censored"} in the prepared data
#'     (if censoring variables are provided).
#' }
#' The non-adherence time is \code{pmin(last_interval, first_deviation, censored_time, na.rm = TRUE)}.
#' The event indicator \code{event_nonadherence} is coded as 0 (censored), 1 (non-adherence),
#' and 2 (competing event / outcome).
#'
#' A stratified cumulative incidence function is fitted using \code{prodlim::prodlim}
#' with \code{Hist(time_nonadherence, event_nonadherence)} and plotted with
#' \code{prodlim::ggprodlim} for cause 1 (non-adherence).
#'
#' @return A \code{ggplot2} object (as returned by \code{prodlim::ggprodlim}) showing the
#'   cumulative incidence of non-adherence (percent) by regime.
#'
#' @examples
#' data(rtmle_object)
#' p <- plot_adherence(rtmle_object)
#' class(p)
#'
#' @seealso \code{\link{intervention_match}}, \code{\link{regime}},
#'   \code{\link{adherence}},
#'   \code{\link{summary_adherence}},
#'   \code{\link{plot_IPW}}, \code{\link[prodlim:prodlim]{prodlim}},
#'   \code{\link[prodlim:ggprodlim]{ggprodlim}},
#'   \code{\link[prodlim:Hist]{Hist}}
#'
#' @importFrom prodlim prodlim ggprodlim Hist
#' @export
plot_adherence <- function(x, regimes = NULL, ...) {
    time_unit <- if (is.null(x$time_unit)) "Time" else x$time_unit
    dt_nonadherence <- adherence(x, regimes = regimes)
    fit_nonadherence <- prodlim::prodlim(Hist(time_nonadherence,event_nonadherence)~regime,
                                         data = dt_nonadherence)
    p <- prodlim::ggprodlim(fit_nonadherence,
                            cause = 1,
                            type = "risk",
                            ylim = c(0,100))
    suppressMessages(p <- p+ggplot2::scale_x_continuous(breaks = x$time_grid, labels = x$time_grid_labels))
    p+ ggplot2::xlab(time_unit)+ ggplot2::ylab("Non-adherence")
}



######################################################################
### plot_adherence.R ends here
