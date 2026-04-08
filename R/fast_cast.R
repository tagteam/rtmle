### fast_cast.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: mar 19 2026 (15:29) 
## Version: 
## Last-Updated: apr  2 2026 (07:54) 
##           By: Thomas Alexander Gerds
##     Update #: 15
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Fast casting from long to wide format on a discrete time grid
#'
#' This function efficiently transforms a long-format \code{data.table}
#' with one row per \code{(id, interval)} into a wide-format table with
#' one row per subject and one column per interval.
#'
#' It is a lightweight alternative to \code{data.table::dcast()} optimized
#' for the specific structure used in discretized time-to-event analyses.
#' @param x A \code{data.table} in long format containing at least the columns
#'   \code{id}, \code{interval}, and a value column.
#'
#' @param name Name of the variable which is transformed from long to wide format.
#' 
#' @param id Character string naming the subject identifier column.
#'
#' @param value_col Character string naming the column in \code{x} that contains
#'   the values to be cast into wide format.
#'
#' @param fill Value used to initialize the output matrix. Cells corresponding
#'   to missing \code{(id, interval)} combinations are filled with this value.
#'
#' @param fun_aggregate Optional function used to aggregate values when multiple
#'   rows exist for the same \code{(id, interval)} combination. If duplicates are
#'   present and \code{fun_aggregate} is \code{NULL}, an error is thrown.
#'
#' @return A \code{data.table} in wide format with:
#'   \describe{
#'     \item{rows}{Unique values of \code{id}}
#'     \item{columns}{One column per interval (named by interval index)}
#'   }
#'
#' @details
#' The function avoids the overhead of \code{dcast()} by directly constructing
#' a matrix and filling it using integer indexing. This makes it substantially
#' faster for large datasets with many intervals.
#'
#' If duplicate \code{(id, interval)} combinations are detected, values are
#' aggregated using \code{fun_aggregate}. Otherwise, values are inserted directly.
#'
#' Column names corresponding to intervals are returned as character strings.
#'
#' Rows are ordered by the first occurrence of each id,
#' and columns are ordered by sorted interval values.
#'
#' @examples
#' library(data.table)
#'
#' x <- data.table(
#'   id = c(1,1,2,2),
#'   interval = c(0,1,0,1),
#'   value = c(10,20,30,40)
#' )
#'
#' fast_cast(x, id = "id")
#'
#' ## With duplicate entries
#' x_dup <- rbind(x, data.table(id = 1, interval = 1, value = 25))
#'
#' fast_cast(x_dup, id = "id", fun_aggregate = mean)
#'
#' @keywords internal
#' @export
fast_cast <- function(x,
                      name = NULL,
                      id,
                      value_col = "value",
                      fill = NA,
                      fun_aggregate = NULL) {
    N <- NULL
    stopifnot(id %in% names(x))
    stopifnot("interval" %in% names(x))
    stopifnot(value_col %in% names(x))

    bycols <- c(id, "interval")
    x <- x[, c(bycols, value_col), with = FALSE]

    # Only aggregate when duplicates are present
    dup <- x[, list(N = .N), by = bycols][N > 1L]
    if (nrow(dup) > 0L) {
        if (is.null(fun_aggregate)) {
            if (length(name)>0){
                stop(paste0("Processing variable ",name,": duplicate (id, interval) combinations found, but fun_aggregate is NULL."))
            }
            stop("Duplicate (id, interval) combinations found, but fun_aggregate is NULL.")
        }
        fun_aggregate <- match.fun(fun_aggregate)
        x <- x[, list(tmp_value = fun_aggregate(.SD[[1L]])),
               by = bycols,
               .SDcols = value_col]
        value_col <- "tmp_value"
    }

    ids <- unique(x[[id]])
    intervals <- sort(unique(x[["interval"]]))

    i <- match(x[[id]], ids)
    j <- match(x[["interval"]], intervals)

    mat <- matrix(fill, nrow = length(ids), ncol = length(intervals))
    ## mat[cbind(i, j)] <- x[[value_col]]
    # can you improve the efficiency of the previous line?
    linear_indices <- (j - 1L) * length(ids) + i
    mat[linear_indices] <- x[[value_col]]
    
    out <- data.table::data.table(tmp_id = ids)
    data.table::setnames(out, "tmp_id", id)
    out <- cbind(out, data.table::as.data.table(mat))
    data.table::setnames(out, c(id, as.character(intervals)))
    out
}

######################################################################
### fast_cast.R ends here
