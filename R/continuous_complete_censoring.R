# Complete the stored censoring table after all baseline and terminal data have
# been added. Keep track of inferred rows so a later preparation can replace
# them if outcome, competing-event, censoring, or baseline data have changed.
continuous_complete_censoring <- function(x) {
    id <- x$names$id
    baseline <- x$data$baseline_data
    if (!is.data.frame(baseline) || !nrow(baseline) || !id %in% names(baseline) ||
        anyNA(baseline[[id]]) || anyDuplicated(baseline[[id]]))
        stop("Add baseline data with nonmissing unique subject identifiers first.")

    censor <- x$long_data$censored_data
    time <- x$names$time
    if (!is.null(censor)) {
        censor <- data.table::copy(data.table::as.data.table(censor))
        if (!time %in% names(censor) && "date" %in% names(censor)) time <- "date"
        if (!all(c(id, time) %in% names(censor)))
            stop("censored_data must contain subject identifiers and a time column.")
        if (!"event" %in% names(censor)) data.table::set(censor, j = "event", value = rep("C", nrow(censor)))
        previous <- x$continuous_inferred_censoring
        # Do not duplicate tauend rows that were inferred earlier
        if (!is.null(previous) && all(c(id, time, "event") %in% names(previous))) {
            index <- match(censor[[id]], previous[[id]])
            inferred <- !is.na(index) & !is.na(censor$event) & censor$event == "tauend" &
                !is.na(censor[[time]]) & censor[[time]] == previous[[time]][index]
            censor <- censor[!inferred]
        }
    } else {
        censor <- data.table::data.table(subject = baseline[[id]][0])
        data.table::setnames(censor, "subject", id)
        data.table::set(censor, j = time, value = numeric())
        data.table::set(censor, j = "event", value = character())
    }

    # Subjects with no recorded terminal event receive an administrative end.
    observed <- censor[[id]]
    for (part in c("outcome_data", "competing_data")) {
        frame <- x$long_data[[part]]
        if (!is.null(frame)) observed <- c(observed, frame[[id]])
    }
    missing <- which(!baseline[[id]] %in% observed)
    inferred <- data.table::data.table(subject = baseline[[id]][missing])
    data.table::setnames(inferred, "subject", id)
    data.table::set(inferred, j = time, value = rep(max(x$time_grid_scale) + 1, length(missing)))
    data.table::set(inferred, j = "event", value = rep("tauend", length(missing)))
    # Create the appropriate columns when there are no censoring records.
    if (!nrow(censor)) censor <- inferred[0]
    x$long_data$censored_data <- data.table::rbindlist(list(censor, inferred), use.names = TRUE, fill = TRUE)
    x$continuous_inferred_censoring <- inferred
    x
}
