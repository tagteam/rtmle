# Keep a diagnostic from the complete event history, before adaptive event
# truncation. It is used only when an analysis requests a static regime.
continuous_treatment_support <- function(data, horizons) {
    id <- A <- NULL
    events <- data.table::copy(data$timevarying_data)
    baseline <- data$baseline_data
    # Compare treatment at each opportunity with the last recorded treatment
    # state, using A_0 before the first event for each person.
    previous <- events[, c(baseline$A_0[match(id[1L], baseline$id)],
                            utils::head(A, -1L)), by = "id"]$V1
    opportunities <- which(events$event == "A" & events$time > 0)
    list(time = events$time[opportunities],
         switched = events$A[opportunities] != previous[opportunities],
         horizons = horizons)
}
