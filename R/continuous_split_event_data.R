# Convert a single table of event marks and current values into the separate
# streams used by add_long_data and the continuous-time preparation pipeline.
continuous_split_event_data <- function(x, data, treatment_events, covariate_events,
                                        time_covariates = NULL) {
    if (!is.data.frame(data) || anyDuplicated(names(data)))
        stop("event_data must be a data frame with unique column names.")
    id <- x$names$id
    time <- if (x$names$time %in% names(data)) x$names$time else
        if ("date" %in% names(data)) "date" else x$names$time
    if (anyDuplicated(c(id, time, "event", "A")) ||
        !all(c(id, time, "event", "A") %in% names(data)))
        stop("event_data must contain distinct subject, time, event, and A columns.")
    continuous_validate_treatment_events(treatment_events)
    if (!is.null(covariate_events) &&
        (!is.character(covariate_events) || !length(covariate_events) ||
         anyNA(covariate_events) || any(!nzchar(covariate_events)) || anyDuplicated(covariate_events)))
        stop("covariate_events must be NULL or distinct, nonempty event labels.")
    terminal_labels <- list(Y = unique(c("Y", x$names$outcome)),
                            D = unique(c("D", x$names$competing)),
                            C = unique(c("C", x$names$censoring)),
                            tauend = "tauend")
    terminals <- unlist(terminal_labels, use.names = FALSE)
    if (anyDuplicated(terminals) || any(c(treatment_events, covariate_events) %in% terminals))
        stop("Terminal event labels must be distinct from treatment and covariate labels.")
    if (is.null(time_covariates)) time_covariates <- setdiff(names(data), c(id, time, "event", "A"))
    if (!is.character(time_covariates) || anyNA(time_covariates) ||
        any(!nzchar(time_covariates)) || anyDuplicated(time_covariates) ||
        any(time_covariates %in% c(id, time, "event", "A")) ||
        !all(time_covariates %in% names(data)))
        stop("time_covariates must name distinct event_data covariate columns other than A.")
    data <- data.table::copy(data.table::as.data.table(data))
    marks <- as.character(data[["event"]])
    if (anyNA(data[[id]]) || anyNA(data[[time]]) || anyNA(marks) ||
        any(!marks %in% c(treatment_events, covariate_events, terminals)))
        stop("event_data contains missing identifiers, times or event labels, or unknown event labels.")
    # A terminal event must stand alone. Nonterminal rows at the same time
    # can describe different parts of one treatment or covariate visit.
    keys <- data[, c(id, time), with = FALSE]
    tied <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
    terminal_tie <- which(tied & marks %in% terminals)
    if (length(terminal_tie)) {
        i <- terminal_tie[1L]
        stop("A terminal event is tied with another event for subject ",
             data[[id]][i], " at time ", data[[time]][i], ".")
    }
    # Covariates may also be measured at a treatment opportunity, even when
    # that event label is not listed separately in covariate_events.
    treatment_rows <- which(marks %in% treatment_events)
    covariate_rows <- which(marks %in% covariate_events)
    covariate_only_rows <- setdiff(covariate_rows, treatment_rows)
    terminal_rows <- which(marks %in% terminals)
    values <- data[["A"]][treatment_rows]
    if (anyNA(values) || !all(as.character(values) %in% c("0", "1")))
        stop("Treatment A must have an observed binary value at each treatment event.")
    observed_A <- data[["A"]][covariate_only_rows]
    if (is.numeric(observed_A) && any(is.nan(observed_A) | is.infinite(observed_A)))
        stop("Treatment A must be finite and binary (0/1) when supplied at a covariate event.")
    observed_A <- observed_A[!is.na(observed_A)]
    if (!all(as.character(observed_A) %in% c("0", "1")))
        stop("Treatment A must be binary (0/1) when supplied at a covariate event.")
    for (variable in time_covariates) {
        values <- data[[variable]]

        if (anyNA(values[covariate_rows]))
            stop("Missing values for ", variable, " at a covariate event.")

        if (is.numeric(values) && any(is.nan(values) | is.infinite(values)))
            stop("Nonfinite values for ", variable, " in event_data.")
    }
    rows <- function(index) data.table::copy(data[index, c(id, time), with = FALSE])
    point <- function(index, variable) {
        out <- rows(index)
        data.table::set(out, j = "value", value = data[[variable]][index])
        # Tied nonterminal rows may describe the same measurement more than
        # once, but two different values for one variable are ambiguous.
        out <- unique(out)
        conflict <- anyDuplicated(out[, c(id, time), with = FALSE])
        if (conflict)
            stop("Conflicting values for ", variable,
                 " for subject ", out[[id]][conflict],
                 " at time ", out[[time]][conflict], ".")
        out
    }
    covariate_points <- function(variable) {
        # Missing values at treatment and terminal events leave the previous
        # measurement in place. Observed values are retained at those times.
        optional_rows <- c(treatment_rows, terminal_rows) ##FIXME: covariate values should be never be considered so checks for it can be ignored
        measured <- optional_rows[!is.na(data[[variable]][optional_rows])]
        point(sort(unique(c(covariate_rows, measured))), variable)
    }
    streams <- c(list(A = point(treatment_rows, "A")),
                 stats::setNames(lapply(time_covariates, covariate_points), time_covariates))
    outcome <- rows(which(marks %in% terminal_labels$Y))
    competing <- rows(which(marks %in% terminal_labels$D))
    censored_rows <- which(marks %in% c(terminal_labels$C, terminal_labels$tauend))
    censored <- rows(censored_rows)
    data.table::set(censored, j = "event", value = ifelse(marks[censored_rows] %in% terminal_labels$C, "C", "tauend"))
    result <- list(outcome_data = outcome, competing_data = competing, censored_data = censored,
                   timevar_data = streams, treatment_events = "A", covariate_events = unique(rows(covariate_rows)))
    # Retain observed A on covariate-only rows for validation once baseline
    # data and date origins are available. These rows are not treatment events.
    recorded <- covariate_only_rows[!is.na(data[["A"]][covariate_only_rows])]
    attr(result, "covariate_event_A") <- point(recorded, "A")
    result
}
