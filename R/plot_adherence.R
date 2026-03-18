### plot_adherence.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: dec 11 2025 (10:23) 
## Version: 
## Last-Updated: mar 14 2026 (08:03) 
##           By: Thomas Alexander Gerds
##     Update #: 25
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Plot cumulative non-adherence by protocol
#'
#' Computes time to first deviation from a treatment regime for each protocol in a
#' \code{speff2trial}-like object \code{x}, allowing for right censoring and competing
#' risks (e.g., death/outcome). The result is plotted as the cumulative incidence
#' (in percent) of non-adherence over follow-up time, stratified by protocol.
#'
#' @param x An object containing protocol-specific adherence information and follow-up
#'   data. Must include at least \code{x$protocols} (a named list where each element has
#'   \code{$intervention_match}), \code{x$followup} (with \code{last_interval}),
#'   \code{x$prepared_data} (optional; used for censoring indicators), and
#'   \code{x$names$censoring}.
#' @param ... Currently unused. Included for future extensions.
#'
#' @details
#' For each protocol, the function restricts to initiators (those with
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
#'   cumulative incidence of non-adherence (percent) by protocol.
#'
#' @examples
#' \dontrun{
#' p <- plot_adherence(x)
#' print(p)
#' }
#'
#' @seealso \code{\link[prodlim:prodlim]{prodlim}}, \code{\link[prodlim:ggprodlim]{ggprodlim}},
#'   \code{\link[prodlim:Hist]{Hist}}
#'
#' @importFrom prodlim prodlim ggprodlim Hist
#' @importFrom data.table data.table :=
#' @export
plot_adherence <- function(x,...){
    time_nonadherence = last_interval = event_nonadherence = value = variable = first_event = event = treatment = NULL
    # find time to first deviation from regime where
    # death is a competing risk and data may be right censored
    dt_nonadherence <- do.call(rbind,lapply(names(x$protocols),function(pro){
        initiators <- x$protocols[[pro]]$intervention_match[,1,drop = TRUE] == 1
        first_deviation <- apply(x$protocols[[pro]]$intervention_match[initiators == 1,,drop = FALSE], 1, function(x) match(0, x))
        if (length(x$names$censoring)>0){
            Cvars <- grep(paste0("^",x$names$censoring,"_[0-9]+$"),names(x$prepared_data),value = TRUE)
            censored_time <- apply(x$prepared_data[initiators == 1,Cvars,with = FALSE], 1, function(x) match("censored", x))
        }
        adherence_data <- cbind(protocol = pro, data.table(x$followup[initiators == 1]),first_deviation = first_deviation,censored_time = censored_time)
        adherence_data[,time_nonadherence := pmin(last_interval,first_deviation,censored_time,na.rm = TRUE)]
        # initialize with 0
        adherence_data[,event_nonadherence := 0]
        # value 2 if not censored (competing risk or outcome)
        adherence_data[is.na(censored_time),event_nonadherence := 2]
        # value 2 when non-adherence is observed
        adherence_data[!is.na(first_deviation),event_nonadherence := 1]
        adherence_data[,data.table::data.table(protocol,time_nonadherence,event_nonadherence)]
    }))
    dt_nonadherence[,protocol := factor(protocol)]
    fit_nonadherence <- prodlim::prodlim(Hist(time_nonadherence,event_nonadherence)~protocol,
                                         data = dt_nonadherence)
    prodlim::ggprodlim(fit_nonadherence,
                       cause = 1,
                       ylim = c(0,100))+ ylab("Non-adherence")
}



######################################################################
### plot_adherence.R ends here
