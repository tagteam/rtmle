### adherence.R ---
#----------------------------------------------------------------------
## Collect non-adherence data by regime
#----------------------------------------------------------------------

resolve_adherence_regimes <- function(x, regimes = NULL) {
    if (is.null(x$regimes) || length(x$regimes) == 0L) {
        stop(
            "rtmle::adherence: No regime has prepared the ",
            "intervention_match table yet.\n",
            "To fix this apply regime() after rtmle_prepare_data() ",
            "or apply run_rtmle() or call intervention_match()."
        )
    }

    regime_names <- names(x$regimes)
    available <- vapply(
        x$regimes,
        function(regime) {
            NROW(regime$intervention_match) > 0L
        },
        logical(1L)
    )
    available_regimes <- regime_names[available]
    if (length(available_regimes) == 0L) {
        stop(
            "rtmle::adherence: No regime has prepared the ",
            "intervention_match table yet.\n",
            "To fix this apply regime() after rtmle_prepare_data() ",
            "or apply run_rtmle() or call intervention_match()."
        )
    }

    if (is.null(regimes)) {
        return(available_regimes)
    }
    if (length(regimes) == 0L ||
        !is.character(regimes) ||
        anyNA(regimes) ||
        any(!nzchar(regimes))) {
        stop("Argument regimes must be a non-empty character vector.")
    }
    regimes <- unique(regimes)
    unavailable <- setdiff(regimes, available_regimes)
    if (length(unavailable) > 0L) {
        prepared_n <- if (is.null(x$prepared_data)) {
            0L
        } else {
            NROW(x$prepared_data)
        }
        if (prepared_n > 0L) {
            stop(
                "The following regimes have no element ",
                "intervention_match yet:\n",
                paste(unavailable, collapse = ", "),
                "\nRun x <- intervention_match(x, regime_name)."
            )
        }
        stop(
            "The object does not contain the prepared data yet.\n",
            "Run x <- prepare_rtmle_data(x)\n",
            "and then x <- intervention_match(x, regime_name)."
        )
    }
    regimes
}

adherence_id_name <- function(x, followup) {
    id_name <- x$names$id
    if (is.null(id_name) ||
        length(id_name) == 0L ||
        isTRUE(is.na(id_name[[1L]]))) {
        id_name <- "id"
    } else {
        id_name <- as.character(id_name[[1L]])
    }
    if (!(id_name %in% names(followup))) {
        stop(
            "The follow-up data must contain the subject identifier column ",
            id_name, "."
        )
    }
    id_name
}

#' Collect non-adherence data by regime
#'
#' Computes the time and event indicator used to describe the first deviation
#' from each treatment regime. Only initiators, namely rows whose first
#' intervention matches the regime, are included.
#'
#' @param x An object containing regime-specific adherence information and
#'   follow-up data. It must include \code{x$regimes} with
#'   \code{$intervention_match} tables and \code{x$followup} with the subject
#'   identifier \code{x$names$id} and a \code{last_interval} column. If
#'   censoring variables are present,
#'   \code{x$prepared_data} and \code{x$names$censoring} are used to identify
#'   censoring times.
#' @param regimes Optional character vector of regime names. If omitted,
#'   all regimes with a prepared \code{intervention_match} table are used.
#' @return A \code{data.table} with one row per regime initiator and the
#'   configured subject identifier column \code{x$names$id}, followed by
#'   \code{regime}, \code{time_nonadherence}, and
#'   \code{event_nonadherence}. The event indicator is 0 for censoring, 1 for
#'   non-adherence, and 2 for a competing event or outcome.
#'   The returned table has class \code{"adherence"}, so
#'   \code{summary(non_adherence)} returns the corresponding
#'   \code{\link{summary_adherence}} table when prepared data are available.
#' @details For each selected regime, the function identifies the first
#'   interval where \code{intervention_match} equals 0 and the first interval
#'   marked \code{"censored"} in the prepared data. The non-adherence time is
#'   \code{pmin(last_interval, first_deviation, censored_time, na.rm = TRUE)}.
#'   This is the data collection step used by \code{\link{plot_adherence}}.
#' @seealso \code{\link{plot_adherence}}, \code{\link{summary_adherence}},
#'   \code{\link{intervention_match}}, \code{\link{regime}}
#' @examples
#' data(rtmle_object)
#' non_adherence <- adherence(rtmle_object, regimes = "Always_A")
#' head(non_adherence)
#' @export
adherence <- function(x, regimes = NULL) {
    regimes <- resolve_adherence_regimes(x, regimes)

    followup <- data.table::as.data.table(x$followup)
    if (!"last_interval" %in% names(followup)) {
        stop("The object must contain follow-up data with a last_interval column.")
    }
    id_name <- adherence_id_name(x, followup)
    empty_adherence <- function() {
        empty <- data.table::data.table(
            regime = character(),
            time_nonadherence = numeric(),
            event_nonadherence = numeric()
        )
        empty[, (id_name) := followup[[id_name]][FALSE]]
        data.table::setcolorder(
            empty,
            c(id_name, "regime", "time_nonadherence", "event_nonadherence")
        )
        empty
    }
    prepared_data <- if (is.null(x$prepared_data)) {
        NULL
    } else {
        data.table::as.data.table(x$prepared_data)
    }
    censoring <- if (is.null(x$names$censoring)) {
        character()
    } else {
        x$names$censoring
    }

    non_adherence <- lapply(regimes, function(regime_name) {
        intervention_match <- x$regimes[[regime_name]]$intervention_match
        if (NCOL(intervention_match) == 0L) {
            return(empty_adherence())
        }
        initiators <- !is.na(intervention_match[, 1L]) &
            intervention_match[, 1L] == 1
        initiator_rows <- which(initiators)
        if (length(initiator_rows) == 0L) {
            return(empty_adherence())
        }

        first_deviation <- apply(
            intervention_match[initiator_rows, , drop = FALSE],
            1L,
            function(values) match(0, values)
        )
        first_deviation <- as.numeric(first_deviation)

        censored_time <- rep(NA_real_, length(initiator_rows))
        if (length(censoring) > 0L && !is.null(prepared_data)) {
            censoring_variables <- unlist(lapply(
                censoring,
                function(variable) {
                    candidates <- names(prepared_data)[startsWith(
                        names(prepared_data),
                        paste0(variable, "_")
                    )]
                    candidates[grepl(
                        "^[0-9]+$",
                        substring(candidates, nchar(variable) + 2L)
                    )]
                }
            ), use.names = FALSE)
            censoring_variables <- unique(censoring_variables)
            if (length(censoring_variables) > 0L) {
                censored_time <- apply(
                    prepared_data[
                        initiator_rows,
                        censoring_variables,
                        with = FALSE
                    ],
                    1L,
                    function(values) {
                        match("censored", as.character(values))
                    }
                )
                censored_time <- as.numeric(censored_time)
            }
        }

        last_interval <- followup[["last_interval"]][initiator_rows]
        time_nonadherence <- pmin(
            last_interval,
            first_deviation,
            censored_time,
            na.rm = TRUE
        )
        event_nonadherence <- numeric(length(initiator_rows))
        event_nonadherence[is.na(censored_time)] <- 2
        event_nonadherence[!is.na(first_deviation)] <- 1

        non_adherence_data <- data.table::data.table(
            regime = regime_name,
            time_nonadherence = time_nonadherence,
            event_nonadherence = event_nonadherence
        )
        non_adherence_data[, (id_name) := followup[[id_name]][initiator_rows]]
        data.table::setcolorder(
            non_adherence_data,
            c(id_name, "regime", "time_nonadherence", "event_nonadherence")
        )
        non_adherence_data
    })

    non_adherence <- data.table::rbindlist(
        non_adherence,
        use.names = TRUE
    )
    non_adherence[, regime := factor(as.character(regime))]
    class(non_adherence) <- c("adherence", class(non_adherence))
    if (!is.null(x$prepared_data)) {
        attr(
            non_adherence,
            "adherence_summary"
        ) <- summary_adherence(x, regimes = regimes)
    }
    non_adherence
}

######################################################################
### adherence.R ends here
