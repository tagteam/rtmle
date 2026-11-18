continuous_invalidate <- function(x) {
    x$continuous_prepared <- x$prepared_data <- x$models <- x$estimate <- x$IC <- NULL
    x$continuous_nuisance <- x$continuous_nuisance_specification <- NULL
    x$continuous_analysis_specification <- NULL
    x$continuous_treatment_support <- NULL
    x$continuous_terminal_times <- x$continuous_baseline_snapshot <- NULL
    x$sequential_outcome_regression <- x$run_time_horizons <- NULL
    x$followup <- NULL
    x$diagnostics <- NULL
    x
}
