### long_to_wide.R ---
#----------------------------------------------------------------------
## Backward-compatible alias for discretize_data()
#----------------------------------------------------------------------

#' Obsolete alias for \code{\link{discretize_data}}
#'
#' \code{long_to_wide()} has been renamed to \code{discretize_data()}.
#' This alias is retained for existing scripts and forwards all arguments to
#' the new function.
#'
#' @param ... Arguments forwarded to \code{\link{discretize_data}}.
#' @return The modified \code{rtmle} object returned by \code{discretize_data}.
#' @seealso \code{\link{discretize_data}}, \code{\link{map_data_to_grid}}
#' @export
long_to_wide <- function(...) {
    .Deprecated("discretize_data", package = "rtmle")
    discretize_data(...)
}

######################################################################
### long_to_wide.R ends here
