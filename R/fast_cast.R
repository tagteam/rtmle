### fast_cast.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: mar 19 2026 (15:29) 
## Version: 
## Last-Updated: mar 19 2026 (15:51) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
fast_cast <- function(x,
                      id,
                      value_col = "value",
                      fill = NA,
                      fun_aggregate = NULL) {
    stopifnot(id %in% names(x))
    stopifnot("interval" %in% names(x))
    stopifnot(value_col %in% names(x))

    bycols <- c(id, "interval")
    x <- x[, c(bycols, value_col), with = FALSE]

    # Only aggregate when duplicates are present
    dup <- x[, list(N = .N), by = bycols][N > 1L]
    if (nrow(dup) > 0L) {
        if (is.null(fun_aggregate)) {
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
    mat[cbind(i, j)] <- x[[value_col]]

    out <- data.table::data.table(tmp_id = ids)
    data.table::setnames(out, "tmp_id", id)
    out <- cbind(out, data.table::as.data.table(mat))
    data.table::setnames(out, c(id, as.character(intervals)))

    out
}

######################################################################
### fast_cast.R ends here
