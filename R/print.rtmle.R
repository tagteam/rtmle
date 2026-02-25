### print.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (08:31) 
## Version: 
## Last-Updated: feb 24 2026 (13:54) 
##           By: Thomas Alexander Gerds
##     Update #: 105
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Printing an object for register data analysis with the targeted minimum loss estimator 
#'
#' This print function informs about the current state of the object and
#' about how to proceed.
#' @param x Object to be printed
#' @param ... Not used for now
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
    cat(sep = "","\nTargeted minimum loss based analysis of longitudinal data on a discretized time scale.\n")
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
    cat(sep = "","\nMinority_Threshold:  ",x$tuning_parameters$minority_threshold)
    cat(sep = "","\nWeight truncation:   ",paste0(as.character(x$tuning_parameters$weight_truncation),sep = ","))
    cat(sep = "","\nPrediction range:    ",paste0(as.character(x$tuning_parameters$prediction_range),sep = ","))
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
        cat(sep = "","\nOutcome data:        n=",NROW(x$data$outcome_data),",p=",format_index(names(x$data$outcome_data)[-1]))
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
    cat(sep = "","\nTODO: Use the function 'prepare_data' to prepare the wide format data.")
    }
    } else {
    cat(sep = "","\nPrepared data:        n=",NROW(x$prepared_data),", p=",(NCOL(prepared_data)-1))
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
    if (length(x$estimate)>0){
        cat("\n\nResults:\n")
        print(x$estimate,...)
    }
    cat("\n")
}

######################################################################
### print.rtmle.R ends here
