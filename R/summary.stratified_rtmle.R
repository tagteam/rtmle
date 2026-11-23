#' Summarize stratified rtmle analyses
#'
#' Combines the estimates and regime contrasts from each stratum.
#'
#' @param object A fitted \code{stratified_rtmle} object.
#' @param ... Arguments passed to \code{\link{summary.rtmle}} for each stratum,
#'   such as \code{targets}, \code{reference}, and \code{digits}.
#' @return A \code{\link[data.table]{data.table}} with a column identifying the
#'   stratum and the estimates from each fit. The column is named after the
#'   stratifying variable, with a numeric suffix if the name is already used.
#' @seealso \code{\link{stratify.rtmle}}, \code{\link{summary.rtmle}}
#' @examples
#' data(rtmle_object)
#' by_sex <- stratify.rtmle(rtmle_object, sex, time_horizon = 1:2)
#' summary(by_sex)
#' @method summary stratified_rtmle
#' @export
summary.stratified_rtmle <- function(object, ...) {
    estimates <- data.table::rbindlist(lapply(object, summary, ...),
                                     use.names = TRUE, fill = TRUE,
                                     idcol = "Stratum")
    data.table::set(estimates, j = "Stratum",
                    value = factor(estimates$Stratum, levels = names(object)))
    variable <- attr(object, "variable", exact = TRUE)
    column_names <- make.unique(c(names(estimates)[-1L], variable))
    data.table::setnames(estimates, 1L, column_names[[length(column_names)]])
    estimates
}

