### clear_rtmle_fit.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Nov 23 2026 (09:36) 
## Version: 
## Last-Updated: Nov 23 2026 (09:36) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
clear_rtmle_fit <- function(x) {
    clear_fitted_components <- function(value) {
        if (!is.list(value)) return(value)
        if (!is.null(names(value))) {
            value <- value[!names(value) %chin% c("fit", "fit_summary")]
        }
        lapply(value, clear_fitted_components)
    }
    x$models <- clear_fitted_components(x$models)
    x$estimate <- NULL
    x$IC <- NULL
    x$sequential_outcome_regression <- NULL
    x$run_time_horizons <- NULL
    x$diagnostics <- NULL
    x$learner <- NULL
    x$unparsed_learner <- NULL
    x$save_fitted_objects <- FALSE
    for (regime_name in names(x$regimes)) {
        x$regimes[[regime_name]]$intervention_probs <- NULL
        x$regimes[[regime_name]]$cumulative_intervention_probs <- NULL
        x$regimes[[regime_name]]$ipw_last_nodes <- NULL
        x$regimes[[regime_name]]$intervention_last_nodes <- NULL
    }
    if (length(x$prepared_data) > 0L &&
        "rtmle_predicted_outcome" %in% names(x$prepared_data)) {
        data.table::set(x$prepared_data,
                        j = "rtmle_predicted_outcome",
                        value = NULL)
    }
    x
}



######################################################################
### clear_rtmle_fit.R ends here
