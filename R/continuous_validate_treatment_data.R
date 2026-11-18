continuous_validate_treatment_data <- function(data) {
    if (is.null(data)) return(invisible(NULL))
    if (all(c("start_date", "end_date") %in% names(data)) ||
        isTRUE(attr(data, "continuous_interval")))
        stop("Treatment A must have one row with an observed value at each treatment opportunity; treatment intervals are not supported.")
    if (!"value" %in% names(data))
        stop("The treatment A table must contain a 'value' column.")
    invisible(NULL)
}
