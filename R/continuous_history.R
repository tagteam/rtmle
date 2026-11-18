continuous_history <- function(data, time_covariates, baseline_covariates,
                                type, lag = NULL, k, exclude_latest_covariate = NULL,
                                model_options = NULL, outcome_name = NULL) {
    previous <- if (k > 1) seq.int(if (is.null(lag)) 1 else max(1, k - lag), k - 1) else integer()
    history <- c(unlist(lapply(c(time_covariates, "time", "event"), function(v) paste0(v, "_", previous))), baseline_covariates)
    if (type == "hazard") history <- setdiff(history, paste0("time_", k - 1))
    if (type == "propensity" && k > 0)
        history <- c(history, paste0(setdiff(c("time", time_covariates), c("A", exclude_latest_covariate)), "_", k))
    history <- intersect(unique(history), names(data))

    # Markov variables use the latest available value, including at baseline.
    # Match the full numeric suffix so a baseline name such as year_2024 is
    # not confused with a different variable's event history.
    indexed_values <- function(variable, columns) {
        prefix <- paste0(variable, "_")
        suffix <- substring(columns, nchar(prefix) + 1L)
        which(startsWith(columns, prefix) & grepl("^[0-9]+$", suffix))
    }
    for (variable in model_options$Markov) {
        indices <- indexed_values(variable, history)
        if (length(indices) > 1L) {
            event_index <- as.numeric(substring(history[indices], nchar(variable) + 2L))
            history <- history[-indices[event_index < max(event_index)]]
        }
    }

    # Rules may restore earlier or globally excluded predictors, but cannot
    # introduce future data or the treatment currently being predicted.
    last_available <- if (type == "propensity") k else max(0L, k - 1L)
    if (type == "baseline") last_available <- 0L
    event_columns <- unlist(lapply(c(time_covariates, "time", "event"), function(variable) {
        paste0(variable, "_", 0:last_available)
    }), use.names = FALSE)
    available <- intersect(unique(c(baseline_covariates, event_columns)), names(data))
    if (type == "propensity")
        available <- setdiff(available, c(paste0("A_", k), paste0(exclude_latest_covariate, "_", k)))
    if (type == "baseline") available <- setdiff(available, "A_0")
    history <- continuous_apply_model_rules(history, available, names(data), model_options,
                                             outcome_name, time_covariates)
    history[vapply(data[, history, with = FALSE], function(v) data.table::uniqueN(v) > 1, logical(1))]
}
