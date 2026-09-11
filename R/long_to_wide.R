### long_to_wide.R ---
#----------------------------------------------------------------------
## Backward-compatible alias for discretize()
#----------------------------------------------------------------------

#' Obsolete alias for \code{\link{discretize}}
#'
#' \code{long_to_wide()} has been renamed to \code{discretize()}.
#' This alias is retained for existing scripts and forwards all arguments to
#' the new function.
#'
#' @param ... Arguments forwarded to \code{\link{discretize}}.
#' @return The modified \code{rtmle} object returned by \code{discretize}.
#' @seealso \code{\link{discretize}}, \code{\link{map_data_to_grid}}
#' @export
long_to_wide <- function(...) {
    .Deprecated("discretize", package = "rtmle")
    discretize(...)
}

######################################################################
### long_to_wide.R ends here
