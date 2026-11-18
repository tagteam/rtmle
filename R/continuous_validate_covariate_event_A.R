# Values of A recorded at covariate events must agree with the treatment
# history; the event itself does not assign treatment.
continuous_validate_covariate_event_A <- function(x) {
    observed <- x$continuous_covariate_event_A
    if (is.null(observed) || !nrow(observed)) return(invisible(NULL))

    id_name <- x$names$id
    time_name <- x$names$time
    observed <- data.table::copy(data.table::as.data.table(observed))
    data.table::setnames(observed, c(id_name, time_name, "value"),
                         c("id", "time", "observed"))
    # A later separate-table call can replace the covariate-event roles.
    # Only validate observations at times that are still covariate events.
    role <- x$continuous_event_roles$covariate
    if (is.null(role)) return(invisible(NULL))
    events <- if (is.character(role)) {
        data.table::rbindlist(lapply(role, function(variable) {
            x$long_data$timevar_data[[variable]][, c(id_name, time_name), with = FALSE]
        }))
    } else data.table::copy(data.table::as.data.table(role))[, c(id_name, time_name), with = FALSE]
    data.table::setnames(events, c(id_name, time_name), c("id", "time"))
    active <- events[observed, on = c("id", "time"), which = TRUE]
    observed <- observed[!is.na(active)]
    if (!nrow(observed)) return(invisible(NULL))

    # Earlier records can precede an unobserved pre-baseline treatment change.
    # At time zero, compare against the initialized baseline value.
    event_time <- observed[["time"]]
    observed <- observed[event_time >= 0]
    if (!nrow(observed)) return(invisible(NULL))

    treatment <- data.table::copy(x$long_data$timevar_data$A)
    data.table::setnames(treatment, c(id_name, time_name, "value"),
                         c("id", "time", "expected"))
    data.table::setorderv(treatment, c("id", "time"))
    values <- treatment[observed, on = c("id", "time"), roll = TRUE]$expected
    if (anyNA(values)) stop("Cannot determine treatment A at a covariate event.")
    mismatch <- which(as.character(observed$observed) != as.character(values))
    if (length(mismatch)) {
        row <- mismatch[1L]
        stop("Treatment A cannot change at a covariate event: subject ",
             observed$id[row], " at time ", observed$time[row], " since baseline.")
    }
    invisible(NULL)
}
