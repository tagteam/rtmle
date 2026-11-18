continuous_check_data <- function(data, time_horizons, time_covariates, baseline_covariates) {
    id <- A <- NULL
    if (!is.numeric(time_horizons) || !length(time_horizons) ||
        any(!is.finite(time_horizons)) || any(time_horizons <= 0) || anyDuplicated(time_horizons))
        stop("time_horizons must contain distinct, finite positive times.")
    for (variables in list(time_covariates, baseline_covariates)) {
        if (!is.character(variables) || !length(variables) || anyNA(variables) ||
            any(!nzchar(variables)) || anyDuplicated(variables))
            stop("Covariate names must be nonempty, unique character vectors.")
    }
    if (!"A" %in% time_covariates || !"A_0" %in% baseline_covariates)
        stop("Include A in time_covariates and A_0 in baseline_covariates.")
    generated_names <- c("time_0", "event_0", "rtmle_pseudo_outcome", "rtmle_treatment_outcome")
    history_names <- grep("_[0-9]+$", baseline_covariates, value = TRUE)
    history_prefix <- sub("_[0-9]+$", "", history_names)
    history_index <- sub("^.*_", "", history_names)
    if (any(baseline_covariates %in% generated_names) ||
        any(history_prefix %in% c("time", "event", "propensity", "survival_censoring", "time_pooled", "event_pooled")) ||
        any(history_prefix %in% time_covariates & history_index != "0"))
        stop("Baseline covariate names must not collide with generated event-history columns.")
    if (!is.list(data) || !all(c("baseline_data", "timevarying_data") %in% names(data)) ||
        !all(vapply(data[c("baseline_data", "timevarying_data")], is.data.frame, logical(1))))
        stop("data must contain baseline_data and timevarying_data data frames.")
    baseline <- data$baseline_data
    events <- data$timevarying_data
    if (!nrow(baseline) || !nrow(events)) stop("Both data frames must contain observations.")
    if (anyDuplicated(names(baseline)) || anyDuplicated(names(events))) stop("Column names must be unique.")
    if (!all(c("id", baseline_covariates) %in% names(baseline))) stop("Missing baseline covariates or id.")
    if (!all(c("id", "time", "event", time_covariates) %in% names(events))) stop("Missing event columns or time-varying covariates.")
    if (anyNA(baseline$id) || anyDuplicated(baseline$id)) stop("Baseline ids must be nonmissing and unique.")
    if (anyNA(events$id) || !setequal(baseline$id, events$id)) stop("Both data frames must contain the same subjects.")
    if (!is.numeric(events$time) || any(!is.finite(events$time)) || any(events$time <= 0))
        stop("Event times must be finite positive numeric values.")
    if (anyDuplicated(as.data.frame(events)[c("id", "time")])) stop("There are ties in event times within a subject.")
    if (anyNA(events$event) || !all(events$event %in% c("A", "L", "C", "Y", "D", "tauend")))
        stop("Invalid event marks; use A, L, C, Y, D, or tauend.")
    check_covariates <- function(frame, variables) {
        if (any(vapply(frame[variables], function(v) anyNA(v) ||
                        (is.numeric(v) && any(!is.finite(v))), logical(1))))
            stop("Covariates must not contain missing (NA) or nonfinite values.")
    }
    check_covariates(as.data.frame(baseline), baseline_covariates)
    check_covariates(as.data.frame(events), time_covariates)
    if (!all(as.character(baseline$A_0) %in% c("0", "1")) ||
        !all(as.character(events$A) %in% c("0", "1"))) stop("Treatment A and A_0 must be binary (0/1).")
    events <- data.table::copy(data.table::as.data.table(events))
    data.table::setorderv(events, c("id", "time"))
    terminal <- !events$event %in% c("A", "L")
    ends <- events[, .I[.N], by = "id"]$V1
    if (sum(terminal) != length(ends) || !all(which(terminal) == ends))
        stop("Each subject must have exactly one terminal event, as their last event.")
    if (any(events$event == "tauend" & events$time < max(time_horizons)))
        stop("End of observation (tauend) must cover the largest time horizon.")
    previous_A <- events[, c(as.character(baseline$A_0[match(id[1], baseline$id)]), utils::head(as.character(A), -1)), by = "id"]$V1
    if (any(as.character(events$A[events$event != "A"]) != as.character(previous_A[events$event != "A"])))
        stop("Treatment may change only at an A event.")
    invisible(NULL)
}
