### plot_IPW.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 26 2026 (09:52) 
## Version: 
## Last-Updated: maj  4 2026 (07:01) 
##           By: Thomas Alexander Gerds
##     Update #: 43
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
#' For each protocol in \code{x$protocols}, extracts cumulative intervention
#' probabilities from \code{cumulative_intervention_probs} and produces
#' boxplots across intervention nodes, restricted to subjects who are:
#' \itemize{
#'   \item \strong{Adherent} at the corresponding intervention decision:
#'   \code{intervention_match[, A_t] == 1}
#'   \item \strong{At risk} at node \code{t}: \code{x$followup$last_interval >= t}
#' }
#'
#' The cumulative probability column is selected from the protocol-specific
#' \code{intervention_last_nodes} index created by \code{\link{run_rtmle}}.
#'
#' @param x An \code{rtmle} object containing:
#'   \itemize{
#'     \item \code{x$protocols}: named list; each protocol has matrices
#'       \code{$cumulative_intervention_probs} and \code{$intervention_match}
#'     \item \code{x$followup}: data frame or data table with columns \code{id} and \code{last_interval}
#'     \item \code{x$intervention_nodes}: integer vector of decision nodes (e.g. \code{c(0,1)})
#'   }
#' @param protocols Character vector of protocol names to include. Default \code{NULL} uses all.
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
                     protocols = NULL
                     ) {
    intervention_node <- intervention_nodes <- used_cumprobs <- NULL
    stopifnot(!is.null(x$protocols),!is.null(x$followup),!is.null(x$intervention_nodes))
    protocol_names <- names(x$protocols)
    if (length(protocol_names) == 0) stop("rtmle::plot_IPW: Object contains no protocols yet. You need to apply 'rtmle::protocol'.") 
    run_protocols <- sapply(protocol_names,function(pn){length(x$protocols[[pn]]$intervention_last_nodes)})
    if (all(run_protocols == 0)) stop("rtmle::plot_IPW: None of the protocols has been fitted to data yet. You need to apply 'rtmle::run_rtmle'.") 
    if (!is.null(protocols)) {
        unavailable_protocols <- setdiff(protocols, protocol_names)
        if (length(unavailable_protocols) > 0) {
            stop(paste0("run_rtmle::plot_IPW: Unavailable protocol(s): ", paste(unavailable_protocols, collapse = ", "), "\nAvailable are: ",paste(names(run_protocols), collapse = ", ")))
        }
        run_protocols <- run_protocols[intersect(names(run_protocols),protocols)]
    }
    plot_dt <- do.call(rbind,lapply(names(run_protocols),function(this_protocol){
        # restrict to those time horizons that have run
        do.call(rbind,lapply(seq_len(run_protocols[[this_protocol]]),function(k){
            outcome_free_and_uncensored <- (x$followup$last_interval >= (k-1))
            if (length(x$names$censoring)>0){
                current_cnode <- as.character(x$prepared_data[[paste0(x$names$censoring,"_",k)]])
                outcome_free_and_uncensored_outcome <- outcome_free_and_uncensored & (current_cnode %in% x$names$uncensored_label)
            }else{
                outcome_free_and_uncensored_outcome <- outcome_free_and_uncensored
            }
            ipos <- x$protocols[[this_protocol]]$intervention_last_nodes[k]
            used_cumprobs <- x$protocols[[this_protocol]]$cumulative_intervention_probs[,ipos]
            if (is.numeric(x$tuning_parameters$weight_truncation)){
                used_cumprobs <- pmax(pmin(used_cumprobs,
                                           x$tuning_parameters$weight_truncation[2]),
                                      x$tuning_parameters$weight_truncation[1])
            }
            intervention_node_name <- paste(x$protocols[[this_protocol]]$intervention_table[time == k-1]$variable,collapse = ",")
            if (nchar(intervention_node_name)>0){
                imatch <- (x$protocols[[this_protocol]]$intervention_match[,intervention_node_name]%in% 1)
            }else{
                imatch <- rep(1,NROW(x$prepared_data))
                imatch[!outcome_free_and_uncensored] <- NA
            }
            subjects_with_weights <- outcome_free_and_uncensored_outcome & as.vector(imatch)
            data.table::data.table(protocol = this_protocol,intervention_node = k - 1,used_cumprobs = used_cumprobs[subjects_with_weights])
        }))
    }))
    plot_dt[, intervention_node := factor(intervention_node, levels = x$intervention_nodes,labels = x$time_grid_labels[x$intervention_nodes+1])]
    missing_values <- plot_dt[,list("missing value" = sum(is.na(used_cumprobs))),by = c("intervention_node","protocol")]
    p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = intervention_node, y = used_cumprobs)) +
        ggplot2::geom_boxplot(outlier.alpha = 0.4) +
        ggplot2::labs(
                     x = "Time",
                     y = "Cumulative intervention probability",
                     title = "Cumulative intervention probabilities among subjects who adhere and are at-risk."
                 ) +
        ggplot2::scale_x_discrete() +
        ggplot2::theme_bw() +
        ggplot2::ylim(c(0,1)) 
        ## ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0))
    if (length(protocol_names)>1) {
        p <- p + ggplot2::facet_grid(. ~ protocol)
    }
    p
}

######################################################################
### plot_IPW.R ends here
