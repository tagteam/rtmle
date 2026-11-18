continuous_validate_treatment_events <- function(events) {
    if (!is.character(events) || !length(events) || anyNA(events) ||
        any(!nzchar(events)) || anyDuplicated(events))
        stop("treatment_events must be a character vector with one or more distinct, nonempty event names.")
    invisible(NULL)
}
