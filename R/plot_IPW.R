### plot_IPW.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 26 2026 (09:52) 
## Version: 
## Last-Updated: apr 25 2026 (09:19) 
##           By: Thomas Alexander Gerds
##     Update #: 39
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
#' @examples
#' x <- list(
#'   protocols = list(Always_A = list(
#'     intervention_last_nodes = c(1, 2),
#'     cumulative_intervention_probs = matrix(c(.9, .8, .7, .6, .9, .7, .5, .3), nrow = 4),
#'     intervention_match = matrix(c(1, 1, 0, 1, 1, 0, 0, 1), nrow = 4,
#'                                 dimnames = list(NULL, c("A_0", "A_1"))),
#'     intervention_table = data.table::data.table(time = c(0, 1), variable = c("A_0", "A_1")))),
#'   followup = data.table::data.table(id = 1:4, last_interval = c(2, 2, 1, 2)),
#'   intervention_nodes = 0:1,
#'   names = list(censoring = "Censored", outcome = "Y",
#'                censored_label = "censored", uncensored_label = "uncensored"),
#'   prepared_data = data.table::data.table(
#'     id = 1:4, Y_1 = 0, Y_2 = 0, Censored_1 = "uncensored",
#'     Censored_2 = c("uncensored", "censored", "uncensored", "uncensored")),
#'   tuning_parameters = list(weight_truncation = c(0, 1)))
#' p <- plot_IPW(x)
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
        do.call(rbind,lapply(x$intervention_nodes[x$intervention_nodes <= run_protocols[[this_protocol]]],function(k){
            outcome_free_and_uncensored <- (x$followup$last_interval >= (k-1))
            if (length(x$names$censoring)>0){
                current_cnode <- as.character(x$prepared_data[[paste0(x$names$censoring,"_",(k+1))]])
                outcome_free_and_uncensored_outcome <- outcome_free_and_uncensored & (current_cnode %in% x$names$uncensored_label)
            }else{
                outcome_free_and_uncensored_outcome <- outcome_free_and_uncensored
            }
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
            Y <- x$prepared_data[[paste0(x$names$outcome,"_",k+1)]]
            subjects_with_weights <- outcome_free_and_uncensored_outcome & as.vector(imatch)
            data.table::data.table(protocol = this_protocol,intervention_node = k,used_cumprobs = used_cumprobs[subjects_with_weights])
        }))
    }))
    plot_dt[, intervention_node := factor(intervention_node, levels = x$intervention_nodes)]
    missing_values <- plot_dt[,list("mising value" = sum(is.na(used_cumprobs))),by = c("intervention_node","protocol")]
    p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = intervention_node, y = used_cumprobs)) +
        ggplot2::geom_boxplot(outlier.alpha = 0.4) +
        ggplot2::labs(
                     x = "Intervention node",
                     y = "Cumulative intervention probability",
                     title = "Cumulative intervention probabilities among intervention-adherent and at-risk subjects"
                 ) +
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
