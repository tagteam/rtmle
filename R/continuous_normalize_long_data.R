# Normalize the existing long-data formats before constructing event histories.
# The caller's raw dates and tables remain available for later preparation.
continuous_normalize_long_data <- function(x) {
    x <- data.table::copy(x)
    baseline <- x$data$baseline_data
    id_name <- x$names$id
    time_name <- x$names$time
    if (!is.data.frame(baseline) || !nrow(baseline) || !id_name %in% names(baseline) ||
        anyNA(baseline[[id_name]]) || anyDuplicated(baseline[[id_name]]))
        stop("Add baseline data with nonmissing unique subject identifiers first.")
    ids <- baseline[[id_name]]
    # Check and describe each input before assembling the event histories.
    describe <- function(frame, label, intervals = TRUE) {
        if (!is.data.frame(frame) || anyDuplicated(names(frame)) || !id_name %in% names(frame))
            stop(label, " must be a data frame with unique column names and subject identifiers.")
        if (anyNA(match(frame[[id_name]], ids)))
            stop("Unknown or missing subject identifiers in ", label, ".")
        interval <- all(c("start_date", "end_date") %in% names(frame))
        if (any(c("start_date", "end_date") %in% names(frame)) && !interval)
            stop(label, " requires both start_date and end_date for interval data.")
        if (interval && !intervals) stop(label, " requires individual event times, not intervals.")
        point <- if (time_name %in% names(frame)) time_name else if ("date" %in% names(frame)) "date" else NULL
        if (!interval && is.null(point)) stop(label, " is missing its time column.")
        columns <- intersect(unique(c(time_name, "date", "start_date", "end_date")), names(frame))
        for (column in columns) {
            values <- frame[[column]]
            if (!is.numeric(values) || is.complex(values) ||
                inherits(values, c("Date", "POSIXt", "difftime")))
                stop(label, " must contain numeric times since baseline; calendar dates are not supported.")
        }
        list(frame = data.table::as.data.table(frame), label = label, interval = interval,
             point = point, columns = columns)
    }
    descriptions <- list()
    for (part in c("outcome_data", "competing_data", "censored_data")) {
        if (!is.null(x$long_data[[part]]))
            descriptions[[part]] <- describe(x$long_data[[part]], part, intervals = FALSE)
    }
    streams <- x$long_data$timevar_data
    if (!is.null(streams) && (!is.list(streams) || is.null(names(streams)) ||
        anyNA(names(streams)) || any(!nzchar(names(streams))) || anyDuplicated(names(streams))))
        stop("timevar_data must be a named list of distinct tables.")
    continuous_validate_treatment_data(streams$A)
    if ("treatment" %in% names(x$continuous_event_roles))
        continuous_validate_treatment_events(x$continuous_event_roles$treatment)
    for (variable in names(streams))
        descriptions[[paste0("stream:", variable)]] <- describe(streams[[variable]], paste0("timevar_data$", variable))
    interval_variables <- names(streams)[vapply(streams, function(frame) {
        all(c("start_date", "end_date") %in% names(frame)) || isTRUE(attr(frame, "continuous_interval"))
    }, logical(1))]
    x$continuous_interval_variables <- interval_variables
    for (role in c("treatment", "covariate")) {
        events <- x$continuous_event_roles[[role]]
        if (is.character(events) && any(events %in% interval_variables))
            stop("Treatment and covariate events must refer to tables with one row per event time, not exposure intervals.")
        if (is.data.frame(events))
            descriptions[[paste0("role:", role)]] <- describe(events, paste0(role, " events"), intervals = FALSE)
    }
    if (!is.null(x$continuous_covariate_event_A) && nrow(x$continuous_covariate_event_A))
        descriptions[["covariate_event_A"]] <- describe(x$continuous_covariate_event_A,
                                                          "Treatment values at covariate events",
                                                          intervals = FALSE)
    check_times <- function(description) {
        frame <- data.table::copy(description$frame)
        for (column in description$columns) {
            values <- as.numeric(frame[[column]])
            if (anyNA(values)) stop("Missing (NA) times in ", description$label, ".")
            infinite_end <- description$interval && identical(column, "end_date")
            if (any(!is.finite(values) & !(infinite_end & values == Inf)))
                stop("Nonfinite times in ", description$label, ".")
            data.table::set(frame, j = column, value = values)
        }
        if (!description$interval && !identical(description$point, time_name))
            data.table::setnames(frame, description$point, time_name)
        frame
    }
    # Validate and copy each table once, before any pre-baseline checks.
    checked <- lapply(descriptions, check_times)

    # Before baseline normalization discards negative times, check whether a
    # covariate-event A value contradicts an earlier observed treatment event.
    # No treatment is assumed before the first recorded treatment event.
    if (!is.null(descriptions[["covariate_event_A"]]) &&
        !is.null(descriptions[["stream:A"]])) {
        observed <- checked[["covariate_event_A"]]
        observed_time <- observed[[time_name]]
        observed <- observed[observed_time < 0]
        if (nrow(observed)) {
            # Resolve the current roles before point-table normalization drops
            # negative times. A later call may have replaced either role.
            role_times <- function(role_name) {
                role <- x$continuous_event_roles[[role_name]]
                keys <- if (is.character(role)) paste0("stream:", role) else paste0("role:", role_name)
                tables <- lapply(keys, function(key) {
                    frame <- checked[[key]]
                    if (is.null(frame)) return(NULL)
                    frame[, c(id_name, time_name), with = FALSE]
                })
                tables <- Filter(Negate(is.null), tables)
                if (!length(tables)) return(observed[0, c(id_name, time_name), with = FALSE])
                unique(data.table::rbindlist(tables))
            }
            treatment <- checked[["stream:A"]]
            treatment_time <- treatment[[time_name]]
            treatment <- treatment[treatment_time < 0]
            treatment_times <- role_times("treatment")
            covariate_times <- role_times("covariate")
            active <- covariate_times[observed, on = c(id_name, time_name), which = TRUE]
            observed <- observed[!is.na(active)]
            active <- treatment_times[treatment, on = c(id_name, time_name), which = TRUE]
            treatment <- treatment[!is.na(active)]
            if (nrow(observed)) {
                treatment_records <- data.table::data.table(
                    subject = treatment[[id_name]], time = treatment[[time_name]],
                    value = as.character(treatment$value), is_treatment = rep(TRUE, nrow(treatment)))
                covariate_records <- data.table::data.table(
                    subject = observed[[id_name]], time = observed[[time_name]],
                    value = as.character(observed$value), is_treatment = rep(FALSE, nrow(observed)))
                records <- data.table::rbindlist(list(treatment_records, covariate_records))
                # A simultaneous treatment assignment precedes its covariate
                # record. Earlier L records can establish the treatment state
                # even when no prior A event was observed.
                data.table::setorderv(records, c("subject", "time", "is_treatment"), c(1, 1, -1))
                previous <- records[, data.table::shift(.SD[[1L]]),
                                    by = "subject", .SDcols = "value"]$V1
                mismatch <- which(!records$is_treatment & !is.na(previous) &
                                  !is.na(records$value) & records$value != previous)
                if (length(mismatch)) {
                    row <- mismatch[1L]
                    stop("Treatment A cannot change at a covariate event: subject ",
                         records$subject[row], " at time ", records$time[row],
                         " since baseline.")
                }
            }
        }
    }

    make_points <- function(subject, times, values) {
        frame <- data.table::data.table(subject = subject, time = times, value = values)
        data.table::setnames(frame, c("subject", "time"), c(id_name, time_name))
        frame
    }
    interval_points <- function(frame, label) {
        if (any(frame$start_date >= frame$end_date)) stop("Each interval end must follow its start in ", label, ".")
        groups <- split(seq_len(nrow(frame)), factor(match(frame[[id_name]], ids), levels = seq_along(ids)))
        pieces <- lapply(seq_along(ids), function(i) {
            rows <- groups[[i]]
            starts <- frame$start_date[rows]
            ends <- frame$end_date[rows]
            ord <- order(starts, ends)
            # Overlapping or touching exposure intervals form a single active
            # period, with value 1 at its start and 0 at its end.
            merged_start <- merged_end <- numeric()
            for (j in ord) {
                n <- length(merged_start)
                if (n && starts[j] <= merged_end[n]) merged_end[n] <- max(merged_end[n], ends[j]) else {
                    merged_start <- c(merged_start, starts[j])
                    merged_end <- c(merged_end, ends[j])
                }
            }
            active <- as.numeric(any(merged_start <= 0 & merged_end > 0))
            positive_starts <- merged_start[merged_start > 0]
            positive_ends <- merged_end[is.finite(merged_end) & merged_end > 0]
            make_points(rep(ids[i], 1L + length(positive_starts) + length(positive_ends)),
                        c(0, positive_starts, positive_ends),
                        c(active, rep(1, length(positive_starts)), rep(0, length(positive_ends))))
        })
        points <- data.table::rbindlist(pieces)
        data.table::setorderv(points, c(id_name, time_name))
        attr(points, "continuous_interval") <- TRUE
        points
    }
    baseline_name <- function(variable) {
        candidates <- c(paste0(variable, "_0"), variable)
        candidates <- candidates[candidates %in% names(baseline)]
        if (length(candidates)) candidates[1L] else NULL
    }
    point_values <- function(frame, variable, label) {
        event_only <- !"value" %in% names(frame)
        if (event_only) data.table::set(frame, j = "value", value = rep(1, nrow(frame)))
        if (anyNA(frame$value) || (is.numeric(frame$value) && any(!is.finite(frame$value))))
            stop("Missing or nonfinite values in ", label, ".")
        keys <- frame[, c(id_name, time_name), with = FALSE]
        duplicate <- duplicated(keys)
        if (any(duplicate)) {
            # Compare values within each subject/time group.
            distinct <- unique(frame[, c(id_name, time_name, "value"), with = FALSE])
            if (anyDuplicated(distinct[, c(id_name, time_name), with = FALSE]))
                stop("Conflicting duplicate timestamps in ", label, ".")
            frame <- frame[!duplicate]
        }
        data.table::setorderv(frame, c(id_name, time_name))
        if (variable == "A") {
            values <- as.character(frame$value)
            if (!all(values %in% c("0", "1")))
                stop("Treatment A must have binary 0/1 values in ", label, ".")
            data.table::set(frame, j = "value", value = as.numeric(values))
            initial_name <- baseline_name("A")
            if (is.null(initial_name)) {
                # Match discrete LOCF at baseline: use the latest previous
                # treatment, with zero when no previous value is available.
                initial <- numeric(length(ids))
                previous <- frame[frame[[time_name]] <= 0]
                if (nrow(previous)) {
                    previous <- previous[!duplicated(previous[[id_name]], fromLast = TRUE)]
                    initial[match(previous[[id_name]], ids)] <- previous$value
                }
            } else {
                initial <- as.character(baseline[[initial_name]])
                if (anyNA(initial) || !all(initial %in% c("0", "1")))
                    stop("Initial treatment A must have binary 0/1 values.")
                initial <- as.numeric(initial)
            }
            # A supplied baseline A or A_0 takes precedence at time zero.
            zero <- make_points(ids, rep(0, length(ids)), initial)
            after <- frame[frame[[time_name]] > 0, c(id_name, time_name, "value"), with = FALSE]
            points <- data.table::rbindlist(list(zero, after))
            data.table::setorderv(points, c(id_name, time_name))
            return(points)
        }
        before <- frame[[time_name]] == 0
        if (event_only) {
            # A supplied baseline value is already the count at time zero.
            # Otherwise, count events at baseline.
            observed_at_baseline <- tabulate(match(frame[[id_name]][before], ids), nbins = length(ids))
            initial_name <- baseline_name(variable)
            if (is.null(initial_name)) {
                initial <- as.numeric(observed_at_baseline)
            } else {
                initial <- as.character(baseline[[initial_name]])
                initial[initial == "TRUE"] <- "1"
                initial[initial == "FALSE"] <- "0"
                initial <- suppressWarnings(as.numeric(initial))
                if (any(!is.finite(initial)) || any(initial < 0) || any(initial != floor(initial)))
                    stop("Event-only baseline counts must be nonnegative whole numbers in ", label, ".")
                if (any(initial < observed_at_baseline))
                    stop("Event-only baseline counts cannot be smaller than the observed baseline events in ", label, ".")
                data.table::set(baseline, j = initial_name, value = initial)
            }
            after <- frame[frame[[time_name]] > 0]
            if (nrow(after)) {
                within_person <- after[, seq_len(.N), by = id_name]$V1
                data.table::set(after, j = "value",
                                value = initial[match(after[[id_name]], ids)] + within_person)
            }
            zero <- make_points(ids, rep(0, length(ids)), initial)
            points <- data.table::rbindlist(list(zero, after[, c(id_name, time_name, "value"), with = FALSE]))
            data.table::setorderv(points, c(id_name, time_name))
            return(points)
        }
        at_zero <- frame[before]
        if (nrow(at_zero)) {
            at_zero <- at_zero[!duplicated(at_zero[[id_name]], fromLast = TRUE)]
            data.table::set(at_zero, j = time_name, value = rep(0, nrow(at_zero)))
        }
        points <- data.table::rbindlist(list(at_zero[, c(id_name, time_name, "value"), with = FALSE],
                                             frame[frame[[time_name]] > 0, c(id_name, time_name, "value"), with = FALSE]))
        data.table::setorderv(points, c(id_name, time_name))
        points
    }
    # Store normalized tables and event roles in the same locations as
    # the original inputs. The downstream merger works only with elapsed times.
    for (name in names(descriptions)) {
        description <- descriptions[[name]]
        frame <- checked[[name]]
        if (startsWith(name, "stream:")) {
            variable <- substring(name, 8L)
            points <- if (description$interval)
                interval_points(frame, description$label) else point_values(frame, variable, description$label)
            if (variable %in% interval_variables) attr(points, "continuous_interval") <- TRUE
            x$long_data$timevar_data[[variable]] <- points
        } else if (startsWith(name, "role:")) {
            if (anyDuplicated(frame[, c(id_name, time_name), with = FALSE]))
                stop("Duplicate event timestamps in ", description$label, ".")
            x$continuous_event_roles[[substring(name, 6L)]] <- frame[frame[[time_name]] >= 0]
        } else if (name == "covariate_event_A") {
            x$continuous_covariate_event_A <- frame
        } else x$long_data[[name]] <- frame
    }
    x
}
