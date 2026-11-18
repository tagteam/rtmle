# One row for each observed interval (start, stop]. Values at its left endpoint
# are held fixed until the next recorded event. Terminal marks are never used.
continuous_censoring_intervals <- function(prepared) {
    events <- data.table::copy(prepared$event_data)
    baseline <- prepared$baseline_data
    if (is.null(events) || is.null(baseline))
        stop("Re-run prepare_rtmle_data before fitting a time-varying censoring model.")
    variables <- prepared$time_covariates
    if (length(intersect(variables, prepared$baseline_covariates)))
        stop("Baseline and current covariates must have distinct names for time-varying censoring.")
    reserved <- c(".rtmle_start", ".rtmle_stop", ".rtmle_censored",
                  ".rtmle_request", ".rtmle_from", ".rtmle_to")
    if (any(c(variables, prepared$baseline_covariates) %in% reserved))
        stop("Covariate names beginning with '.rtmle_' are reserved for censoring intervals.")
    initial <- paste0(variables, "_0")
    if (!all(initial %in% names(baseline)))
        stop("Time-varying censoring requires baseline values named: ", paste(initial, collapse = ", "), ".")
    data.table::setorderv(events, c("id", "time"))
    first <- !duplicated(events$id)
    starts <- c(0, utils::head(events$time, -1L))
    starts[first] <- 0
    intervals <- events[, "id", with = FALSE]
    data.table::set(intervals, j = ".rtmle_start", value = starts)
    data.table::set(intervals, j = ".rtmle_stop", value = events$time)
    data.table::set(intervals, j = ".rtmle_censored", value = as.integer(events$event == "C"))
    for (variable in variables) {
        values <- events[, data.table::shift(.SD[[1L]]), by = "id", .SDcols = variable]$V1
        initial <- baseline[[paste0(variable, "_0")]][match(events$id[first], baseline$id)]
        categorical <- is.factor(values) || is.factor(initial)
        if (categorical) {
            values <- as.character(values)
            initial <- as.character(initial)
        }
        values[first] <- initial
        if (categorical) values <- factor(values, levels = unique(values))
        data.table::set(intervals, j = variable, value = values)
    }
    baseline[intervals, on = "id"]
}
