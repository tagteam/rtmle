#' Add event-time data to a continuous-time rtmle analysis
#'
#' @param x An rtmle object initialized with \code{continuous_time = TRUE}.
#' @inheritParams add_long_data
#' @param outcome_data A data frame with subject identifiers and outcome times.
#'   Times are durations since baseline; only subjects with an outcome
#'   need a row. Durations need to be called \code{time}.
#' @param censored_data A data frame with subject identifiers and censoring
#'   times, as numeric durations since baseline. An optional \code{event} column
#'   distinguishes censoring from the end of follow-up. Durations need to be called \code{time}.
#' @param competing_data A data frame with subject identifiers and competing
#'   event times, as durations since baseline. Only subjects with a
#'   competing event need a row. Durations need to be called \code{time}.
#' @param timevar_data A named list of treatment and covariate tables, each
#'   with subject identifiers and times since baseline. The treatment
#'   table also needs a \code{value} column. Durations need to be called \code{time}.
#' @param event_data An optional wide event table with identifier and time
#'   columns named as configured in \code{rtmle_init}, an \code{event} column,
#'   an \code{A} column, and any time-varying covariate columns. Event labels
#'   identify treatment, covariate measurement and terminal events. 
#' @param time_covariates Names of covariate columns in \code{event_data},
#'   excluding \code{A}. With \code{NULL}, infer them from the other columns;
#'   \code{character(0)} selects none.
#' @param treatment_events A nonempty character vector naming distinct tables in
#'   \code{timevar_data}, each with one row per treatment event time, or event
#'   labels in \code{event_data}. See examples.
#' @param covariate_events Covariate-measurement times: names of tables in
#'   \code{timevar_data} with one row per measurement time, a separate
#'   identifier/time data frame, or
#'   \code{NULL} for none. With \code{event_data}, use event labels or
#'   \code{NULL}. Treatment and measurement events may coincide.
#' @details 
#' Use numeric times since baseline. Add baseline data with
#' \code{\link{add_baseline_data}}, then use \code{\link{prepare_rtmle_data}}
#' to assemble the event histories.
#' 
#' Record binary treatment \code{A} at every treatment opportunity. For
#' baseline \code{A}, a supplied baseline value takes precedence. Otherwise,
#' RTMLE uses the latest treatment at or before time zero, or 0 if none is
#' recorded.
#'
#' The time-varying covariate data is supplied either in \code{event_data}
#' which has the data fully prepared or in \code{timevar_data}
#' which prepares the data for the user.
#' With \code{event_data}, specify both \code{treatment_events} and
#' \code{covariate_events}. Use \code{NULL} if there are no separate covariate
#' events. Supply all selected covariates on rows marked as covariate events.
#' At other events, missing covariates keep their previous values. Covariates
#' measured at a treatment visit are treated as available before assignment.
#'
#' In separate tables for \code{event_data}, a covariate table without a \code{value} column
#' represents a cumulative event count. Exposure intervals give the current
#' binary state. Declare measurement times separately.
#' 
#' Use \code{"Y"}, \code{"D"}, \code{"C"}, and \code{"tauend"} for outcome,
#' competing event, censoring, and administrative end of follow-up. The outcome,
#' competing-event, and censoring names set in \code{rtmle_init()} also work.
#' A terminal event must have a distinct time for that subject. If none is
#' recorded, preparation adds \code{"tauend"} one time unit after the last
#' time in the grid.
#' 
#' @return The rtmle object, ready for \code{prepare_rtmle_data}.
#' @seealso \code{\link{run_rtmle}}, \code{\link{prepare_continuous_data}}
#' @examples
#' x <- rtmle_init(c(0, 5), name_id = "id", name_time = "time",
#'                 name_outcome = "Y", name_competing = "D",
#'                 name_censoring = "C", continuous_time = TRUE)
#' x <- add_baseline_data(x, data.frame(id = 1:4, A = c(0, 1, 0, 1),
#'                                      L1 = c(10, 20, 30, 40),
#'                                      L2 = c(100, 200, 300, 400)))
#'
#' # One table: event labels identify treatment and covariate updates.
#' events <- data.frame(
#'     id = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4),
#'     time = c(1, 2, 4, 1, 2, 4, 1, 2, 4, 1, 2),
#'     event = c("A", "visit", "Y", "A", "visit", "D", "A", "visit", "C", "A", "visit"),
#'     A = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1),
#'     L1 = c(11, 12, 13, 21, 22, 22, 31, 32, 32, 41, 42),
#'     L2 = c(101, 102, 103, 201, 202, 202, 301, 302, 302, 401, 402)
#' )
#' x_events <- add_continuous_data(x, event_data = events,
#'                                 treatment_events = "A",
#'                                 covariate_events = "visit",
#'                                 time_covariates = c("L1", "L2"))
#' x_events <- prepare_rtmle_data(x_events, min_events = 0)
#'
#' # Equivalent separate tables: each point table has an id, time and value.
#' outcomes <- data.frame(id = 1, time = 4)
#' competing <- data.frame(id = 2, time = 4)
#' censoring <- data.frame(id = 3, time = 4)
#' treatment <- data.frame(id = 1:4, time = 1, value = c(1, 1, 0, 1))
#' # Both covariates change at treatment visits and other visits; subject 1
#' # also has measurements recorded at the time of the outcome.
#' measurement_id <- c(1, 1, 1, 2, 2, 3, 3, 4, 4)
#' measurement_time <- c(1, 2, 4, 1, 2, 1, 2, 1, 2)
#' l1 <- data.frame(id = measurement_id, time = measurement_time,
#'                  value = c(11, 12, 13, 21, 22, 31, 32, 41, 42))
#' l2 <- data.frame(id = measurement_id, time = measurement_time,
#'                  value = c(101, 102, 103, 201, 202, 301, 302, 401, 402))
#' # Terminal-time measurements are not separate covariate events.
#' covariate_visits <- data.frame(id = 1:4, time = 2)
#' x_tables <- add_continuous_data(x, outcome_data = outcomes,
#'                                 censored_data = censoring,
#'                                 competing_data = competing,
#'                                 timevar_data = list(A = treatment, L1 = l1, L2 = l2),
#'                                 treatment_events = "A",
#'                                 covariate_events = covariate_visits)
#' x_tables <- prepare_rtmle_data(x_tables, min_events = 0)
#' # Subject 4 has no terminal event; preparation infers tauend at time 6.
#' x_tables$long_data$censored_data
#' all.equal(x_events$prepared_data, x_tables$prepared_data)
#' @export
add_continuous_data <- function(x, outcome_data, censored_data, competing_data,
                                timevar_data, treatment_events, covariate_events,
                                event_data, time_covariates = NULL, ...) {
    if (!inherits(x, "rtmle") || !isTRUE(x$continuous_time))
        stop("Initialize x with rtmle_init(continuous_time = TRUE).")
    if (!missing(event_data)) {
        if (!missing(outcome_data) || !missing(censored_data) || !missing(competing_data) ||
            !missing(timevar_data) || length(list(...)))
            stop("event_data cannot be combined with separate long-data inputs in the same call.")
        if (missing(treatment_events) || missing(covariate_events))
            stop("Declare both treatment_events and covariate_events for event_data; use NULL for no covariate events.")
        # Convert the combined format once, then use the same path as separate tables.
        event_streams <- continuous_split_event_data(x, event_data, treatment_events,
                                                     covariate_events, time_covariates)
        recorded_A <- attr(event_streams, "covariate_event_A", exact = TRUE)
        x$long_data$timevar_data <- NULL
        x$continuous_covariate_event_A <- NULL
        x <- do.call(add_continuous_data, c(list(x = x), event_streams))
        x$continuous_covariate_event_A <- recorded_A
        return(x)
    }
    if (!missing(time_covariates)) stop("time_covariates is only used with event_data.")
    if (!missing(treatment_events)) continuous_validate_treatment_events(treatment_events)
    args <- list(x = x, ...)
    if (!missing(outcome_data)) args$outcome_data <- outcome_data
    if (!missing(censored_data)) args$censored_data <- censored_data
    if (!missing(competing_data)) args$competing_data <- competing_data
    if (!missing(timevar_data)) args$timevar_data <- timevar_data
    x <- do.call(add_long_data, args)
    if (!missing(censored_data) && !is.null(censored_data)) x$continuous_inferred_censoring <- NULL
    continuous_validate_treatment_data(x$long_data$timevar_data$A)
    for (part in intersect(names(args), c("outcome_data", "censored_data", "competing_data")))
        if (is.data.frame(args[[part]]) && nrow(args[[part]]) == 0L)
            x$long_data[[part]] <- NULL
    roles <- x$continuous_event_roles
    if (is.null(roles)) roles <- list()
    if (!missing(treatment_events)) roles["treatment"] <- list(data.table::copy(treatment_events))
    if (!missing(covariate_events)) roles["covariate"] <- list(data.table::copy(covariate_events))
    x$continuous_event_roles <- roles
    continuous_invalidate(x)
}
