### add_baseline_data.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (08:54)
## Version:
## Last-Updated: feb 23 2026 (13:36) 
##           By: Thomas Alexander Gerds
##     Update #: 19
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
### Code:
#' Add baseline data to an rtmle object
#' 
#' Adds a baseline data set to an existing \code{rtmle} object.
#'
#' @param x An object of class \code{"rtmle"}, typically created by
#'   \code{\link{rtmle_init}}.
#' @param data A data frame with one row per subject and a subject identifier
#'   column matching \code{x$names$id}. If follow-up starts at a subject-specific
#'   date, include that start-of-follow-up column here.
#' @param ... Not used.
#' @return The modified \code{rtmle} object.
#' @seealso \code{\link{add_long_data}}, \code{\link{add_wide_data}},
#'   \code{\link{long_to_wide}}, \code{\link{prepare_rtmle_data}}
#' @examples
#' x <- rtmle_init(time_grid = 0:2, name_id = "id", name_outcome = "Y")
#' baseline <- data.frame(id = 1:3,
#'                        age = c(55, 60, 65),
#'                        sex = c(0, 1, 0))
#' x <- add_baseline_data(x, baseline)
#' x$data$baseline_data
#' @export
add_baseline_data <- function(x,data,...){
    if (!is.data.frame(data) || !(x$name$id %in% names(data)))
        stop(paste0("'data' needs to be a data.frame (or data.table or tibble) with a variable called ",x$names$id))
    if (!(x$names$id %in% names(data)))
        stop(paste0("'data' needs to contain the subject identifier variable: '",x$names$id,"'"))
    data.table::setDT(data)
    if (any(duplicated(data[[x$names$id]]))) stop("Duplicated values of subject id variable ",x$names$id," are not allowed in baseline data.")
    x$data$baseline_data <- data.table::copy(data)
    x
}

######################################################################
### add_baseline_data.R ends here
