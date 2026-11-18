# Merge the same separate data streams accepted by add_long_data, retaining
# their observed times rather than mapping them onto discrete intervals.
continuous_data_from_long <- function(x) {
    id <- time <- event <- priority <- value <- NULL
    baseline <- x$data$baseline_data
    if (!is.data.frame(baseline) || !nrow(baseline))
        stop("Add baseline data with add_baseline_data before preparing continuous data.")
    baseline <- data.table::copy(data.table::as.data.table(baseline))
    id_name <- x$names$id
    if (!id_name %in% names(baseline) || anyNA(baseline[[id_name]]) ||
        anyDuplicated(baseline[[id_name]])) stop("Baseline subject identifiers must be nonmissing and unique.")
    if (anyDuplicated(names(baseline))) stop("Baseline column names must be unique.")
    if (id_name != "id" && "id" %in% names(baseline))
        stop("Baseline column id conflicts with the canonical subject identifier.")
    data.table::setnames(baseline, id_name, "id")
    streams <- x$long_data$timevar_data
    if (is.null(streams)) streams <- list()
    continuous_validate_treatment_data(streams$A)
    if (is.list(streams) && !is.data.frame(streams) && !"A" %in% names(streams) &&
        any(c("A", "A_0") %in% names(baseline))) {
        empty_A <- data.table::data.table(subject = baseline$id[0], time = numeric(), value = numeric())
        data.table::setnames(empty_A, c("subject", "time"), c(id_name, x$names$time))
        streams$A <- empty_A
    }
    variables <- names(streams)
    if (!is.list(streams) || is.null(variables) || anyNA(variables) ||
        any(!nzchar(variables)) || anyDuplicated(variables) || !"A" %in% variables)
        stop("timevar_data must be a named list of distinct covariates including A.")
    if (any(variables %in% c("id", "time", "event")))
        stop("Covariate names id, time, and event are reserved for continuous histories.")
    horizon <- max(x$time_grid_scale)
    event_only <- vapply(streams, function(frame) is.data.frame(frame) && !"value" %in% names(frame), logical(1))

    # Normalize names on copies, including the customary date-column fallback.
    normalize <- function(frame, label, terminal = FALSE) {
        if (!is.data.frame(frame)) stop(label, " must be a data frame.")
        frame <- data.table::copy(data.table::as.data.table(frame))
        if (anyDuplicated(names(frame))) stop(label, " has duplicated column names.")
        if (!id_name %in% names(frame)) stop(label, " is missing subject identifiers.")
        if (all(c("start_date", "end_date") %in% names(frame)))
            stop(label, ": continuous histories require rows with an id, time, and value, not exposure intervals.")
        time_name <- if (x$names$time %in% names(frame)) x$names$time else
            if ("date" %in% names(frame)) "date" else NULL
        if (is.null(time_name)) stop(label, " is missing its time column.")
        times <- frame[[time_name]]
        if (!is.numeric(times) || inherits(times, c("Date", "POSIXt", "difftime")) ||
            any(!is.finite(times)) || any(times < 0) || (terminal && any(times == 0)))
            stop(label, ": times must be finite numeric durations since baseline",
                 if (terminal) " and strictly positive." else " and nonnegative.")
        matched <- match(frame[[id_name]], baseline$id)
        if (anyNA(matched)) stop("Unknown subject identifiers in ", label, ".")
        if (anyDuplicated(frame[, c(id_name, time_name), with = FALSE]))
            stop("Duplicated times within a subject in ", label, ".")
        if (terminal && anyDuplicated(frame[[id_name]]))
            stop(label, " must contain at most one terminal time per subject.")
        data.table::set(frame, j = id_name, value = baseline$id[matched])
        data.table::setnames(frame, c(id_name, time_name), c("id", "time"))
        frame
    }
    updates <- lapply(variables, function(variable) {
        frame <- normalize(streams[[variable]], paste0("timevar_data$", variable))
        if (!"value" %in% names(frame)) data.table::set(frame, j = "value", value = rep(1, nrow(frame)))
        if (anyNA(frame$value) || (is.numeric(frame$value) && any(!is.finite(frame$value))))
            stop("Missing or nonfinite update values in ", variable, ".")
        frame[, list(id, time, value)]
    })
    names(updates) <- variables

    # Reconcile baseline columns with observations at time zero and give each
    # categorical history one set of levels for subsequent carry-forward.
    for (variable in variables) {
        initial_name <- paste0(variable, "_0")
        if (variable %in% names(baseline)) {
            if (initial_name %in% names(baseline)) {
                if (!identical(as.character(baseline[[variable]]), as.character(baseline[[initial_name]])))
                    stop("Baseline columns ", variable, " and ", initial_name, " conflict.")
                data.table::set(baseline, j = variable, value = NULL)
            } else data.table::setnames(baseline, variable, initial_name)
        }
        at_zero <- updates[[variable]][time == 0]
        zero_match <- match(baseline$id, at_zero$id)
        zero_value <- at_zero$value[zero_match]
        if (!initial_name %in% names(baseline)) {
            if (event_only[[variable]]) zero_value[is.na(zero_match)] <- 0
            else if (anyNA(zero_match)) stop("Missing initial values for ", variable, "; supply ", initial_name, " at baseline or updates at time zero.")
            data.table::set(baseline, j = initial_name, value = zero_value)
        } else if (nrow(at_zero)) {
            supplied <- !is.na(zero_match)
            if (anyNA(baseline[[initial_name]]) ||
                any(as.character(baseline[[initial_name]][supplied]) != as.character(zero_value[supplied])))
                stop("Time-zero updates conflict with supplied baseline values for ", variable, ".")
        }
        initial <- baseline[[initial_name]]
        if (anyNA(initial) || (is.numeric(initial) && any(!is.finite(initial))))
            stop("Missing or nonfinite initial values for ", variable, ".")
        values <- updates[[variable]]$value
        if (variable == "A") {
            if (!all(c(as.character(initial), as.character(values)) %in% c("0", "1")))
                stop("Treatment A and its initial values must be binary (0/1).")
            data.table::set(baseline, j = initial_name, value = as.numeric(as.character(initial)))
            data.table::set(updates[[variable]], j = "value", value = as.numeric(as.character(values)))
        } else if (is.factor(initial) || is.character(initial) || is.factor(values) || is.character(values)) {
            # Use a shared level set so carry-forward also preserves new categories.
            levels <- unique(c(levels(initial), levels(values), as.character(initial), as.character(values)))
            data.table::set(baseline, j = initial_name, value = factor(as.character(initial), levels = levels))
            data.table::set(updates[[variable]], j = "value", value = factor(as.character(values), levels = levels))
        } else if (!((is.numeric(initial) || is.logical(initial)) &&
                     (is.numeric(values) || is.logical(values)))) {
            stop("Initial and update values for ", variable, " must have compatible numeric or categorical types.")
        }
    }

    # Keep the first terminal event per subject after rejecting tied events.
    terminals <- list()
    for (part in c("outcome_data", "competing_data", "censored_data")) {
        frame <- x$long_data[[part]]
        if (is.null(frame)) next
        frame <- normalize(frame, part, terminal = TRUE)
        if (!nrow(frame)) next
        if (part == "censored_data") {
            if ("event" %in% names(frame)) {
                marks <- as.character(frame$event)
                marks[marks %in% x$names$censoring] <- "C"
                if (anyNA(marks) || !all(marks %in% c("C", "tauend")))
                    stop("censored_data event marks must be C or tauend.")
            } else marks <- rep("C", nrow(frame))
        } else {
            mark <- if (part == "outcome_data") "Y" else "D"
            configured <- if (part == "outcome_data") x$names$outcome else x$names$competing
            if ("event" %in% names(frame) &&
                (anyNA(frame$event) || !all(as.character(frame$event) %in% c(mark, configured))))
                stop(part, " contains inconsistent terminal event marks.")
            marks <- rep(mark, nrow(frame))
        }
        terminals[[part]] <- data.table::data.table(id = frame$id, time = frame$time, event = marks)
    }
    terminal <- data.table::rbindlist(terminals)
    terminal_records <- data.table::copy(terminal)
    terminal_tie_error <- function(record) {
        stop("Event times tied for subject ", as.character(record$id[1L]),
             " at time ", format(record$time[1L]),
             ": a terminal event cannot coincide with another event.", call. = FALSE)
    }
    if (nrow(terminal)) {
        tied <- duplicated(terminal[, c("id", "time"), with = FALSE])
        if (any(tied)) terminal_tie_error(terminal[which(tied)[1L]])
        data.table::setorderv(terminal, c("id", "time"))
        terminal <- terminal[!duplicated(id), list(id, time, event)]
    }
    missing <- baseline$id[!baseline$id %in% terminal$id]
    inferred <- data.table::data.table(id = missing, time = rep(horizon + 1, length(missing)),
                                      event = rep("tauend", length(missing)))
    terminal <- data.table::rbindlist(list(terminal, inferred), use.names = TRUE)
    terminal_records <- data.table::rbindlist(list(terminal_records, inferred), use.names = TRUE)

    roles <- x$continuous_event_roles
    if (!is.list(roles) || !all(c("treatment", "covariate") %in% names(roles)))
        stop("Declare treatment_events and covariate_events explicitly; use NULL for no covariate events.")
    continuous_validate_treatment_events(roles$treatment)
    interval_variables <- variables[vapply(streams, function(frame) isTRUE(attr(frame, "continuous_interval")), logical(1))]
    interval_variables <- union(interval_variables, x$continuous_interval_variables)
    role_times <- function(role, label) {
        if (is.null(role)) return(data.table::data.table(id = baseline$id[0], time = numeric()))
        if (is.character(role)) {
            if (!length(role) || anyNA(role) || anyDuplicated(role) || !all(role %in% variables))
                stop(label, " must name distinct tables with one row per event time.")
            if (any(role %in% interval_variables))
                stop(label, " cannot use interval boundaries as opportunities; supply explicit event times.")
            times <- data.table::rbindlist(lapply(updates[role], function(frame) frame[, list(id, time)]))
        } else {
            times <- normalize(role, label)[, list(id, time)]
        }
        unique(times[time > 0])
    }
    treatment_times <- role_times(roles$treatment, "treatment_events")
    covariate_times <- role_times(roles$covariate, "covariate_events")
    event_times <- unique(data.table::rbindlist(list(treatment_times, covariate_times)))
    tied <- terminal_records[event_times, on = c("id", "time"), nomatch = 0L]
    if (nrow(tied)) terminal_tie_error(tied)
    # An opportunity must include the observed treatment even if it did not
    # change. Treatment measurements at or after a terminal event are ignored.
    relevant_times <- treatment_times[time < terminal$time[match(id, terminal$id)]]
    if (nrow(relevant_times) &&
        anyNA(updates$A[relevant_times, on = c("id", "time"), which = TRUE]))
        stop("Treatment A must have an observed value at every treatment opportunity.")
    # Covariates can also be measured at a terminal event without creating
    # an additional covariate-event row in the merged history.
    measurement_times <- unique(data.table::rbindlist(list(
        treatment_times, covariate_times, terminal[, list(id, time)])))
    histories <- list()
    for (variable in variables) {
        points <- updates[[variable]][time > 0]
        endpoints <- terminal$time[match(points$id, terminal$id)]
        # The terminal event precedes any simultaneous treatment decision.
        points <- if (variable == "A") points[time < endpoints] else points[time <= endpoints]
        if (!variable %in% interval_variables && nrow(points) &&
            anyNA(measurement_times[points, on = c("id", "time"), which = TRUE]))
            stop("Undeclared measurement events in ", variable,
                 "; every observation must occur at a treatment, covariate, or terminal event time.")
        initial <- data.table::data.table(id = baseline$id, time = 0,
                                          value = baseline[[paste0(variable, "_0")]])
        history <- data.table::rbindlist(list(initial, points))
        data.table::setorderv(history, c("id", "time"))
        if (variable == "A") {
            previous <- history[, data.table::shift(as.character(value)), by = "id"]$V1
            changes <- which(!is.na(previous) & as.character(history$value) != previous)
            if (length(changes) &&
                anyNA(treatment_times[history[changes], on = c("id", "time"), which = TRUE]))
                stop("Treatment A may change only at declared treatment opportunity times.")
        }
        histories[[variable]] <- history
    }
    marked_times <- function(times, mark) {
        data.table::data.table(id = times$id, time = times$time, event = rep(mark, nrow(times)))
    }
    # Merge simultaneous treatment and covariate roles into one A event.
    # Carry all covariates to each event, including terminal events.
    schedule <- data.table::rbindlist(list(terminal, marked_times(treatment_times, "A"),
                                          marked_times(covariate_times, "L")))
    schedule[, priority := match(event, c("Y", "D", "C", "tauend", "A", "L"))]
    data.table::setorderv(schedule, c("id", "time", "priority"))
    schedule <- schedule[!duplicated(schedule[, list(id, time)])]
    end_times <- terminal$time[match(schedule$id, terminal$id)]
    schedule <- schedule[time <= end_times, list(id, time, event)]
    for (variable in variables) {
        carried <- histories[[variable]][schedule, on = c("id", "time"), roll = TRUE]$value
        data.table::set(schedule, j = variable, value = carried)
    }
    list(data = list(baseline_data = baseline, timevarying_data = schedule),
         time_covariates = variables,
         baseline_covariates = setdiff(names(baseline), "id"))
}
