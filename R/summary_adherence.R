### summary_adherence.R ---
#----------------------------------------------------------------------
## Summarize adherence and end-of-follow-up events by regime
#----------------------------------------------------------------------

#' Summarize adherence and end-of-follow-up events by regime
#'
#' For each selected regime and intervention interval, counts the people who
#' initiated that regime, still adhere to it, and have their follow-up end
#' because of the outcome, a competing event (death), or censoring.
#'
#' @param x An \code{rtmle} object containing prepared data, follow-up
#'   information, and regime-specific \code{intervention_match} tables.
#' @param regimes Optional character vector of regime names. If omitted,
#'   all regimes with a prepared \code{intervention_match} table are used.
#' @return A \code{data.table} with one row per selected regime and interval.
#'   The columns are \code{regime}, \code{time_node},
#'   \code{n_initiated}, \code{n_adherent}, \code{n_outcome},
#'   \code{n_death}, and \code{n_censored}. The first two counts refer to
#'   regime initiators; the last three counts identify the cause of the end
#'   of follow-up in that interval.
#' @details
#' A person is counted as adherent at an interval when the cumulative
#' \code{intervention_match} value is 1. End-of-follow-up events are counted
#' only in the interval in which they occur. The outcome, competing-event,
#' and censoring variables use suffixes corresponding to the interval's time
#' node plus one, as in \code{x$prepared_data}.
#' @seealso \code{\link{adherence}}, \code{\link{intervention_match}},
#'   \code{\link{plot_adherence}}
#' @examples
#' data(rtmle_object)
#' adherence_summary <- summary_adherence(rtmle_object)
#' adherence_summary
#' @export
summary_adherence <- function(x, regimes = NULL) {
    regimes <- resolve_adherence_regimes(x, regimes)

    if (is.null(x$prepared_data)) {
        stop(
            "summary_adherence requires x$prepared_data to count ",
            "outcome, death, and censoring events."
        )
    }
    prepared_data <- data.table::as.data.table(x$prepared_data)
    followup <- data.table::as.data.table(x$followup)
    if (!"last_interval" %in% names(followup)) {
        stop("The object must contain follow-up data with a last_interval column.")
    }

    outcome <- x$names$outcome
    competing <- x$names$competing
    censoring <- x$names$censoring
    censored_label <- x$names$censored_label
    if (is.null(censored_label) || length(censored_label) == 0L) {
        censored_label <- "censored"
    }
    initiator_rows_for <- function(intervention_match) {
        if (NCOL(intervention_match) == 0L) {
            return(integer())
        }
        which(!is.na(intervention_match[, 1L]) &
              intervention_match[, 1L] == 1)
    }
    count_event <- function(prefix, event_node, labels, initiator_rows) {
        if (is.null(prefix) ||
            length(prefix) == 0L ||
            is.null(labels) ||
            length(labels) == 0L) {
            return(0L)
        }
        variable <- paste0(prefix, "_", event_node)
        if (!(variable %in% names(prepared_data))) {
            return(0L)
        }
        event <- as.character(prepared_data[[variable]]) %in%
            as.character(labels)
        if (length(initiator_rows) == 0L) {
            return(0L)
        }
        at_end <- followup[["last_interval"]][initiator_rows] ==
            event_node - 1L
        as.integer(sum(event[initiator_rows] & at_end, na.rm = TRUE))
    }

    summaries <- lapply(regimes, function(regime_name) {
        intervention_match <- x$regimes[[regime_name]]$intervention_match
        n_nodes <- NCOL(intervention_match)
        if (n_nodes == 0L) {
            return(data.table::data.table(
                regime = character(),
                time_node = integer(),
                n_initiated = integer(),
                n_adherent = integer(),
                n_outcome = integer(),
                n_death = integer(),
                n_censored = integer()
            ))
        }

        initiator_rows <- initiator_rows_for(intervention_match)
        time_node <- seq_len(n_nodes) - 1L
        n_initiated <- length(initiator_rows)
        n_adherent <- if (length(initiator_rows) == 0L) {
            rep.int(0L, n_nodes)
        } else {
            colSums(
                intervention_match[initiator_rows, , drop = FALSE] == 1,
                na.rm = TRUE
            )
        }
        event_node <- time_node + 1L
        n_outcome <- vapply(
            event_node,
            count_event,
            integer(1L),
            prefix = outcome,
            labels = "1",
            initiator_rows = initiator_rows
        )
        n_death <- vapply(
            event_node,
            count_event,
            integer(1L),
            prefix = competing,
            labels = "1",
            initiator_rows = initiator_rows
        )
        n_censored <- vapply(
            event_node,
            count_event,
            integer(1L),
            prefix = censoring,
            labels = censored_label,
            initiator_rows = initiator_rows
        )

        data.table::data.table(
            regime = regime_name,
            time_node = time_node,
            n_initiated = rep.int(n_initiated, n_nodes),
            n_adherent = as.integer(n_adherent),
            n_outcome = n_outcome,
            n_death = n_death,
            n_censored = n_censored
        )
    })

    data.table::rbindlist(summaries, use.names = TRUE)
}

#' @rdname summary_adherence
#' @param object An object returned by \code{\link{adherence}}.
#' @param ... Not used.
#' @export
#' @method summary adherence
summary.adherence <- function(object, ...) {
    adherence_summary <- attr(object, "adherence_summary", exact = TRUE)
    if (is.null(adherence_summary)) {
        stop(
            "This adherence object has no summary. ",
            "Create it from an object with prepared data using adherence(x)."
        )
    }
    adherence_summary
}

######################################################################
### summary_adherence.R ends here
