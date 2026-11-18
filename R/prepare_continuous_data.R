#' Prepare event histories for continuous-time estimation
#'
#' Low-level preparation shared by the continuous-time rtmle workflow.
#' Input tables are copied and ordered by subject and time.
#' @param data A list with canonical \code{baseline_data} and
#'   \code{timevarying_data} frames, with columns \code{id}, \code{time},
#'   \code{event}, \code{A}, and covariates. This is an internal event-history
#'   representation; the analysis workflow accepts separate long tables.
#' @param time_covariates Names of time-varying covariates, including \code{A}.
#' @param baseline_covariates Names of baseline predictors, including \code{A_0}.
#' @param min_events Minimum number of nonterminal events needed to retain an index.
#' @param last_non_terminal_event Optional nonnegative event-index cutoff.
#' @param time_horizons Distinct positive analysis horizons.
#' @param verbose Whether to report adaptive event cutoffs.
#' @return A list with wide event histories, censoring data, covariate names,
#'   and horizon-specific pooling information.
#' @keywords internal
prepare_continuous_data <- function(data, time_horizons, time_covariates,
                                    baseline_covariates,
                                    min_events = 40, last_non_terminal_event = NULL,
                                    verbose = FALSE) {
    event_number <- event <- N <- last_event <- is_censored <- id <- NULL
    continuous_check_data(data, time_horizons, time_covariates, baseline_covariates)
    if (!is.numeric(min_events) || length(min_events) != 1L || !is.finite(min_events) || min_events < 0)
        stop("min_events must be a nonnegative number.")
    if (!is.null(last_non_terminal_event) &&
        (!is.numeric(last_non_terminal_event) || length(last_non_terminal_event) != 1L ||
         !is.finite(last_non_terminal_event) || last_non_terminal_event < 0 ||
         last_non_terminal_event != floor(last_non_terminal_event)))
        stop("last_non_terminal_event must be a nonnegative integer.")
    baseline <- data.table::copy(data.table::as.data.table(data$baseline_data))
    events <- data.table::copy(data.table::as.data.table(data$timevarying_data))
    data.table::set(baseline, j = "A_0", value = as.numeric(as.character(baseline$A_0)))
    data.table::set(events, j = "A", value = as.numeric(as.character(events$A)))
    data.table::set(events, j = "event", value = factor(events$event, levels = c("A", "L", "C", "Y", "D", "tauend")))
    data.table::setorderv(events, c("id", "time"))
    # Keep all recorded updates for censoring models, before outcome-history pooling.
    full_events <- data.table::copy(events)
    events[, event_number := seq_len(.N), by = "id"]
    # Choose the last event index supported by enough subjects at each horizon.
    # Retain terminal events even when later nonterminal events are omitted.
    cutoff <- vapply(time_horizons, function(th) {
        if (!is.null(last_non_terminal_event)) return(as.numeric(last_non_terminal_event))
        counts <- events[time < th & event %in% c("A", "L"), .N, by = "event_number"]
        eligible <- counts[N > min_events, event_number]
        if (length(eligible)) max(eligible) else 0
    }, numeric(1))
    keep_through <- cutoff[which.max(time_horizons)]
    events <- events[event_number <= keep_through | !event %in% c("A", "L")]
    events[, event_number := seq_len(.N), by = "id"]
    info <- data.table::data.table(time_horizon = time_horizons,
                                  last_non_terminal_event = cutoff)
    info[, last_event := vapply(seq_along(time_horizons), function(j) {
        outcome_indices <- events[time <= time_horizons[j] & event == "Y", event_number]
        max(1, min(cutoff[j] + 1, if (length(outcome_indices)) max(outcome_indices) else 1))
    }, numeric(1))]
    info[, is_censored := vapply(time_horizons, function(th) any(events$event == "C" & events$time <= th), logical(1))]
    if (verbose) message("Last event indices: ", paste(info$last_event, collapse = ", "))
    # Sequential regression uses one column per event index. If a terminal
    # event occurs beyond a horizon's last index, keep it in pooled columns.
    wide <- data.table::dcast(events, id ~ event_number, value.var = c("time", "event", time_covariates))
    wide <- wide[baseline, on = "id"]
    data.table::set(wide, j = "event_0", value = rep("A", nrow(wide)))
    data.table::set(wide, j = "time_0", value = rep(0, nrow(wide)))
    for (k in unique(info$last_event)) {
        pooled <- events[!event %in% c("A", "L") & event_number >= k, list(id, time, event)]
        data.table::setnames(pooled, c("time", "event"), paste0(c("time_pooled_", "event_pooled_"), k))
        wide <- pooled[wide, on = "id"]
    }
    censoring_data <- if (any(info$is_censored))
        events[!event %in% c("A", "L")][baseline, on = "id"] else NULL
    structure(list(wide_data = wide, info = info, time_covariates = time_covariates,
                   baseline_covariates = baseline_covariates,
                   event_data = full_events, baseline_data = baseline,
                   data_marginal_censoring = censoring_data), class = "prepare_data_continuous")
}
