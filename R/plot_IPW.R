### plot_IPW.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 26 2026 (09:52) 
## Version: 
## Last-Updated: feb 26 2026 (10:37) 
##           By: Thomas Alexander Gerds
##     Update #: 20
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
#' For each protocol in `x$protocols`, this function extracts cumulative intervention
#' probabilities from `cumulative_intervention_probs` and produces boxplots across
#' intervention nodes, restricted to subjects who are:
#' \itemize{
#'   \item \strong{Adherent} at the corresponding intervention decision: `intervention_match[, A_t] == 1`
#'   \item \strong{At risk} at node `t`: `x$followup$last_interval >= t`
#' }
#'
#' \strong{Default column mapping (requested):}
#' \itemize{
#'   \item Node 0 uses cumulative probability column \code{"Censored_1"} (paired with \code{"A_0"})
#'   \item Node 1 uses cumulative probability column \code{"Censored_2"} (paired with \code{"A_1"})
#' }
#'
#' You can optionally plot \emph{all} columns of `cumulative_intervention_probs`
#' (e.g., including \code{A_0}, \code{A_1}, etc.) instead of only the censored columns.
#'
#' @param x An object containing:
#'   \itemize{
#'     \item \code{x$protocols}: named list; each protocol has matrices
#'       \code{$cumulative_intervention_probs} and \code{$intervention_match}
#'     \item \code{x$followup}: data.frame/data.table with columns \code{id} and \code{last_interval}
#'     \item \code{x$intervention_nodes}: integer vector of decision nodes (e.g. \code{c(0,1)})
#'   }
#' @param protocols Character vector of protocol names to include. Default \code{NULL} uses all.
#' @return A \code{ggplot} object.
#' @export
plot_IPW <- function(
                     x,
                     protocols = NULL
                     ) {
    intervention_nodes <- used_cumprobs <- NULL
    stopifnot(!is.null(x$protocols),!is.null(x$followup),!is.null(x$intervention_nodes))
    protocol_names <- names(x$protocols)
    if (!is.null(protocols)) {
        missing_protocols <- setdiff(protocols, protocol_names)
        if (length(missing_protocols) > 0) stop("Unknown protocols: ", paste(missing_protocols, collapse = ", "))
        protocol_names <- protocols
    }else{
        protocols <- protocol_names
    }
    plot_dt <- do.call(rbind,lapply(protocol_names,function(this_protocol){
        do.call(rbind,lapply(x$intervention_nodes,function(k){
            outcome_free_and_uncensored <- (x$followup$last_interval >= (k-1))
            ipos <- x$protocols[[this_protocol]]$intervention_last_nodes[k+1]
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
            subjects_with_weights <- outcome_free_and_uncensored & as.vector(imatch)
            data.table::data.table(protocol = this_protocol,intervention_node = k,used_cumprobs = used_cumprobs[subjects_with_weights])
        }))
    }))
    plot_dt[, intervention_node := factor(intervention_node, levels = x$intervention_nodes)]
    p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = intervention_node, y = used_cumprobs)) +
        ggplot2::geom_boxplot(outlier.alpha = 0.4) +
        ggplot2::labs(
                     x = "Intervention node",
                     y = "Cumulative intervention probability",
                     title = "Cumulative intervention probabilities among intervention-adherent and at-risk subjects"
                 ) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    if (length(protocol_names)>1) {
        p <- p + ggplot2::facet_grid(. ~ protocol)
    }
    p
}

######################################################################
### plot_IPW.R ends here
