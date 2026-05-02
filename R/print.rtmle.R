### print.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (08:31) 
## Version: 
## Last-Updated: apr 10 2026 (15:51) 
##           By: Thomas Alexander Gerds
##     Update #: 189
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Print an rtmle object
#'
#' Prints the current state of an \code{rtmle} object and indicates the next
#' analysis steps that are still missing.
#'
#' @param x Object to print.
#' @param ... Not used.
#' @return Invisibly returns \code{NULL}.
#' @seealso \code{\link{rtmle_init}}, \code{\link{prepare_rtmle_data}},
#'   \code{\link{run_rtmle}}, \code{\link{summary.rtmle}}
#' @examples
#' x <- rtmle_init(time_grid = 0:2,
#'                 name_id = "id",
#'                 name_outcome = "Y")
#' print(x)
#' @method print rtmle
#' @export
print.rtmle <- function(x, ...) {
    format_index <- function(x, max_show = 5) {
        n <- length(x)
        if (n <= max_show) {
            paste(x, collapse = ", ")
        } else {
            paste0(x[1], ", ", x[2], ", ..., ", x[n])
        }
    }
    cat(sep = "","\nTargeted minimum loss-based analysis of longitudinal data on a discretized time scale.\n")
    # Package version
    if (!is.null(x$version)) {
        cat(sep = "","\nPackage version:     ", as.character(x$version))
    }
    cat(sep = "","\nSubject id variable: ", x$names$id)
    cat(sep = "","\nOutcome variable:    ", x$names$outcome)
    if (length(x$names$censoring)>0){
        cat(sep = "","\nCensoring variable:  ", x$names$censoring)
    }
    if (length(x$names$competing)>0){
        cat(sep = "","\nCompeting risk:      ", x$names$competing)
    }
    cat(sep = "","\nIntervention nodes:  ", format_index(x$intervention_nodes))
    # Tuning parameters
    cat(sep = "","\nMinority threshold:  ",x$tuning_parameters$minority_threshold)
    cat(sep = "","\nWeight truncation:   (",paste0(as.character(x$tuning_parameters$weight_truncation),collapse = ","),")")
    cat(sep = "","\nPrediction range:    (",paste0(as.character(x$tuning_parameters$prediction_range),collapse = ","),")")
    # Data
    if (length(x$data$baseline_data) > 0){
    cat(sep = "","\nBaseline data:       n=",NROW(x$data$baseline_data),", p=",(NCOL(x$data$baseline_data)-1))
    }else{
        cat(sep = "","\nTODO: The object contains no baseline_data yet. Add them with the function rtmle::add_baseline_data")
    }
    if (length(x$intervention_nodes)>1){
        if (length(x$data$timevar_data) > 0){
    cat(sep = "","\nTimevar data:        ",format_index(names(x$data$timevar_data)))
        }else{
            cat(sep = "","\nTODO: The object contains no time varying data yet. Add them with rtmle::add_wide_data in discretized form or with rtmle::add_long_data followed by rtmle::long_to_wide.")
        }
    }
    if (length(x$data$outcome_data) > 0){
    cat(sep = "","\nOutcome data:        n=",NROW(x$data$outcome_data),", p=",format_index(names(x$data$outcome_data)[-1]))
    }else{
        cat(sep = "","\nTODO: The object contains no outcome data yet. Add them with rtmle::add_wide_data in discretized form or with rtmle::add_long_data followed by rtmle::long_to_wide.")
    }
    # Target trial protocols
    if (length(x$protocols) > 0) {
    cat(sep = "","\nProtocols:           ", paste(names(x$protocols), collapse = ", "))
    } else {
        cat(sep = "","\nTODO: The object contains no protocols. Add them with the function 'protocol'.")
    }
    
    if (length(x$targets) > 0) {
        for (t in names(x$targets)) {
            cat(sep = "","\nTarget:              ", t, " [",paste(x$targets[[t]]$protocols, collapse = ", "),"]")
        }
        if (length(x$prepared_data) == 0) {
            cat(sep = "","\nTODO: Use the function 'prepare_rtmle_data' to prepare the wide format data.")
        }
    } else {
        cat(sep = "","\nPrepared data:       n=",NROW(x$prepared_data),", p=",(NCOL(x$prepared_data)-1))
    }
    if (length(x$models) == 0) {
        cat(sep = "","\nTODO: Use the function 'model_formula' to initialize the formula for the nuisance parameter models.")
    } else {
        count_formulas <- function(obj) {
            if (is.list(obj)) {
                sum(vapply(obj, count_formulas, integer(1)))
            } else {
                as.integer(is.character(obj) && length(obj) == 1)
            }
        }
        cat(sep = "","\nNumber of models:    ", count_formulas(x$models))
    }
    summarize_learners <- function(x,called_from_summarize_learners = FALSE){
        if (is.null(x$learners)){
            # --------------------------------
            # SINGLE LEARNER
            # --------------------------------
            learner_args <- setdiff(names(x),c("name","fun"))
            if (length(x[[learner_args]]) == 0){
                args_string <- ""
            } else{
                args <- x[[learner_args]]
                args_string <- paste(
                    paste0(
                        names(args), "=",
                        vapply(args, function(a){
                            paste0(deparse(a), collapse = "")
                        }, character(1))
                    ),
                    collapse = ", "
                )
            }
            if (called_from_summarize_learners){
                if (is.null(args_string)){
                    return(paste0("                       +", x$fun))
                } else {
                    return(paste0("                       +", x$fun, " (", args_string, ")"))
                }
            }else{
                if (args_string[[1]] == "" || length(args_string) == 0){
                    x$fun
                }else{
                    paste0("", x$fun, " (", args_string, ")")
                }
            }
        } else {
            # --------------------------------
            # SUPERLEARNER
            # --------------------------------
            folds <- x$args$folds
            header <- paste0("SuperLearner (", folds, "-fold CV)")
            learner_lines <- unlist(
                lapply(x$learners, function(l){summarize_learners(l,called_from_summarize_learners = TRUE)})
            )
            paste0(header,"\n", paste0(learner_lines,collapse = "\n"))
        }
    }
    if (length(x$learner)>0){
        cat(sep = "","\nLearner:             ",summarize_learners(x = x$learner))
    }else{
        cat(sep = "","\nTODO: Use the function 'run_rtmle' to estimate the nuisance parameter models and the target parameter.")
    }
    if(length(x$run_time_horizons)>0){
        cat(sep = "","\nFitted time horizons:   ", format_index(x$run_time_horizons))
    }
    if (length(x$diagnostics)>0){
        cat("\nDiagnostics:")
        for (diag in names(x$diagnostics)){
            cat("\n                    ",paste0(diag," (n=",NROW(x$diagnostics[[diag]]),")"))
        }
    }else{
        cat("\nDiagnostics:         No warnings (yet)")
    }
    if (length(x$estimate)>0){
        cat("\n\nResults:\n")
        sx = summary(x)[,c("Target","Protocol","Target_parameter","Time_horizon","Estimator","Estimate (CI_95)"),with = FALSE]
        print(sx)
    }
    cat("\n")
}

######################################################################
### print.rtmle.R ends here
