continuous_prepare_rtmle <- function(x, ...) {
    if (!all(c("treatment", "covariate") %in% names(x$continuous_event_roles)))
        stop("Declare treatment_events and covariate_events in add_continuous_data (use NULL for no covariate events).")
    continuous_validate_treatment_events(x$continuous_event_roles$treatment)
    continuous_validate_treatment_data(x$long_data$timevar_data$A)
    supplied <- list(...)
    removed <- intersect(names(supplied), c("start_followup_date", "baseline_lookback"))
    if (length(removed))
        stop("Continuous-time preparation no longer supports: ", paste(removed, collapse = ", "),
             ". Supply times as numeric durations since baseline.")
    settings <- x$continuous_preparation_settings
    settings[c("start_followup_date", "baseline_lookback")] <- NULL
    settings <- utils::modifyList(if (is.null(settings)) list() else settings, supplied)
    if ("marginal_censoring" %in% names(settings))
        stop("Specify marginal_censoring in censoring_learner when running the analysis, not when preparing data.")
    # Complete missing administrative ends, then assemble and widen the event histories.
    # ADD TAUEND TO MISSING CENSORING DATA (infer it)
    x <- continuous_complete_censoring(x)
    normalized <- continuous_normalize_long_data(x)
    input <- continuous_data_from_long(normalized)
    continuous_validate_covariate_event_A(normalized)
    support <- continuous_treatment_support(input$data, x$time_grid_scale[-1])
    prepared <- do.call(prepare_continuous_data,
                        c(input, list(time_horizons = x$time_grid_scale[-1]), settings))
    x <- continuous_invalidate(x)
    x$continuous_preparation_settings <- settings
    x$continuous_prepared <- prepared
    x$continuous_treatment_support <- support
    x$continuous_baseline_snapshot <- data.table::copy(x$data$baseline_data)
    x$prepared_data <- data.table::copy(x$continuous_prepared$wide_data)
    data.table::setnames(x$prepared_data, "id", x$names$id)
    x$names$name_time_covariates <- x$continuous_prepared$time_covariates
    x$names$name_baseline_covariates <- x$continuous_prepared$baseline_covariates
    events <- input$data$timevarying_data
    terminal <- !events$event %in% c("A", "L")
    x$continuous_terminal_times <- data.table::data.table(id = events$id[terminal], time = events$time[terminal])
    # Existing rtmle output uses interval indices, so preserve actual event
    # times separately and express follow-up in those indices.
    last_interval <- pmin(length(x$time_grid) - 2L,
                          findInterval(events$time[terminal], x$time_grid_scale, left.open = TRUE) - 1L)
    x$followup <- data.table::data.table(id = events$id[terminal], last_interval = last_interval)
    data.table::setnames(x$followup, "id", x$names$id)
    data.table::setnames(x$continuous_terminal_times, "id", x$names$id)
    x
}
