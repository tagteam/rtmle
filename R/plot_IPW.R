### plot_IPW.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 26 2026 (09:52) 
## Version: 
## Last-Updated: maj 28 2026 (14:58) 
##           By: Thomas Alexander Gerds
##     Update #: 48
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Boxplots of cumulative intervention probabilities among adherent, at-risk subjects
#'
#' For each regime in \code{x$regimes}, extracts cumulative intervention
#' probabilities from \code{cumulative_intervention_probs} and produces
#' boxplots across intervention nodes, restricted to subjects who are:
#' \itemize{
#'   \item \strong{Adherent} at the corresponding intervention decision:
#'   \code{intervention_match[, A_t] == 1}
#'   \item \strong{At risk} at node \code{t}: \code{x$followup$last_interval >= t}
#' }
#'
#' The cumulative probability column is selected from the regime-specific
#' \code{ipw_last_nodes} index created by \code{\link{run_rtmle}}.
#'
#' @param x An \code{rtmle} object containing:
#'   \itemize{
#'     \item \code{x$regimes}: named list; each regime has matrices
#'       \code{$cumulative_intervention_probs} and \code{$intervention_match}
#'     \item \code{x$followup}: data frame or data table with columns \code{id} and \code{last_interval}
#'     \item \code{x$intervention_nodes}: integer vector of decision nodes (e.g. \code{c(0,1)})
#'   }
#' @param regimes Character vector of regime names to include. Default \code{NULL} uses all.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{run_rtmle}}, \code{\link{plot_adherence}},
#'   \code{\link{plot.rtmle}}
#' @examples
#' data(rtmle_object)
#' p <- plot_IPW(rtmle_object)
#' class(p)
#' @export
plot_IPW <- function(
                     x,
                     regimes = NULL
                     ) {
    time_node <- intervention_nodes <- used_cumprobs <- NULL
    time_unit <- if (is.null(x$time_unit)) "Time" else x$time_unit
    stopifnot(!is.null(x$regimes),!is.null(x$followup),!is.null(x$intervention_nodes))
    regime_names <- names(x$regimes)
    if (length(regime_names) == 0) stop("rtmle::plot_IPW: Object contains no regimes yet. You need to apply 'rtmle::regime'.") 
    run_regimes <- sapply(regime_names,function(pn){length(x$regimes[[pn]]$ipw_last_nodes)})
    if (all(run_regimes == 0)) stop("rtmle::plot_IPW: None of the regimes has been fitted to data yet. You need to apply 'rtmle::run_rtmle'.") 
    if (!is.null(regimes)) {
        unavailable_regimes <- setdiff(regimes, regime_names)
        if (length(unavailable_regimes) > 0) {
            stop(paste0("run_rtmle::plot_IPW: Unavailable regime(s): ", paste(unavailable_regimes, collapse = ", "), "\nAvailable are: ",paste(names(run_regimes), collapse = ", ")))
        }
        run_regimes <- run_regimes[intersect(names(run_regimes),regimes)]
    }
    plot_dt <- do.call(rbind,lapply(names(run_regimes),function(this_regime){
        # restrict to those time horizons that have run
        do.call(rbind,lapply(seq_len(run_regimes[[this_regime]]),function(k){
            outcome_free_and_uncensored <- (x$followup$last_interval >= (k-1))
            if (length(x$names$censoring)>0){
                current_cnode <- as.character(x$prepared_data[[paste0(x$names$censoring,"_",k)]])
                outcome_free_and_uncensored_outcome <- outcome_free_and_uncensored & (current_cnode %in% x$names$uncensored_label)
            }else{
                outcome_free_and_uncensored_outcome <- outcome_free_and_uncensored
            }
            ipos <- x$regimes[[this_regime]]$ipw_last_nodes[k]
            used_cumprobs <- x$regimes[[this_regime]]$cumulative_intervention_probs[,ipos]
            if (is.numeric(x$tuning_parameters$weight_truncation)){
                used_cumprobs <- pmax(pmin(used_cumprobs,
                                           x$tuning_parameters$weight_truncation[2]),
                                      x$tuning_parameters$weight_truncation[1])
            }
            intervention_node_name <- x$regimes[[this_regime]]$intervention_last_nodes[[paste0("node_",k-1)]]
            if (!is.na(intervention_node_name)){
                imatch <- (x$regimes[[this_regime]]$intervention_match[,intervention_node_name]%in% 1)
            }else{
                imatch <- rep(1,NROW(x$prepared_data))
                imatch[!outcome_free_and_uncensored] <- NA
            }
            subjects_with_weights <- outcome_free_and_uncensored_outcome & as.vector(imatch)
            data.table::data.table(
                            regime = this_regime,
                            time_node = k - 1,
                            used_cumprobs = used_cumprobs[subjects_with_weights]
                        )
        }))
    }))
    plot_dt[, time_node := factor(time_node, levels = x$intervention_nodes,labels = x$time_grid_labels[x$intervention_nodes+1])]
    missing_values <- plot_dt[,list("missing value" = sum(is.na(used_cumprobs))),by = c("time_node","regime")]
    p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = time_node, y = used_cumprobs)) +
        ggplot2::geom_boxplot(outlier.alpha = 0.4) +
        ggplot2::labs(
                     x = time_unit,
                     y = "Cumulative intervention probability",
                     title = "Cumulative intervention probabilities among subjects who adhere and are at-risk."
                 ) +
        ggplot2::scale_x_discrete() +
        ggplot2::theme_bw() +
        ggplot2::ylim(c(0,1)) 
        ## ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0))
    if (length(regime_names)>1) {
        p <- p + ggplot2::facet_grid(. ~ regime)
    }
    p
}

######################################################################
### plot_IPW.R ends here
