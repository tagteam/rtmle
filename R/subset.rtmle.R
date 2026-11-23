#' Select participants in a prepared rtmle object
#'
#' Restricts the data to selected participants and clears fitted results.
#' Use \code{\link{run_rtmle}} to fit the analysis on the selected data.
#'
#' @param x A discrete-time \code{rtmle} object with data prepared by
#'   \code{\link{prepare_rtmle_data}}.
#' @param subset Logical expression evaluated in \code{x$prepared_data}.
#'   Rows where the expression is \code{TRUE} are selected.
#' @param ids Participant IDs to select. Supply either \code{subset} or
#'   \code{ids}.
#' @param ... Unused.
#' @return An unfitted \code{rtmle} object with the selected participants,
#'   model formulas, targets and regimes.
#' @seealso \code{\link{run_rtmle}}, \code{\link{cheap_bootstrap}}
#' @examples
#' \dontrun{
#' data(rtmle_object)
#' women <- subset(rtmle_object, sex == 0)
#' women <- run_rtmle(women, time_horizon = 2)
#' summary(women)
#' }
#' @method subset rtmle
#' @export
subset.rtmle <- function(x, subset, ids, ...) {
    if (isTRUE(x$continuous_time)) {
        stop("subset() is only available for discrete-time analyses.")
    }
    id_name <- x$names$id
    data <- x$prepared_data
    if (NROW(data) == 0L || !id_name %in% names(data)) {
        stop("Prepare the data with prepare_rtmle_data() before calling subset().")
    }
    if (missing(subset) == missing(ids)) {
        stop("Supply either a logical subset expression or participant ids.")
    }
    if (missing(ids)) {
        keep <- eval(substitute(subset), data, parent.frame())
        if (!is.logical(keep) || !length(keep) %in% c(1L, nrow(data))) {
            stop("The subset expression must return a logical vector ",
                 "of length one or one value per participant.")
        }
        rows <- which(rep_len(keep, nrow(data)))
    } else {
        rows <- which(data[[id_name]] %in% ids)
    }
    if (length(rows) == 0L) {
        stop("No participants were selected.")
    }
    x <- subset_rtmle_rows(x, rows)
    ids <- x$prepared_data[[id_name]]
    select_ids <- function(data) {
        if (is.data.frame(data)) {
            keep <- data[[id_name]] %in% ids
            return(data[keep, , drop = FALSE])
        }
        if (is.list(data)) return(lapply(data, select_ids))
        data
    }
    x$data <- select_ids(x$data)
    x$long_data <- select_ids(x$long_data)
    if (!is.null(x[["progress"]]$discretized_grid)) {
        x$progress$discretized_grid <- select_ids(x$progress$discretized_grid)
    }
    x
}

