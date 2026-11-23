#' Plot stratified rtmle analyses
#'
#' Plots fitted estimates with one row per stratum and one column per target.
#'
#' @param x A fitted \code{stratified_rtmle} object.
#' @param ... Arguments passed to \code{autoplot.stratified_rtmle}, including
#'   \code{targets}, \code{target_parameter}, \code{reference}, \code{conf_int},
#'   \code{xlim}, and \code{ylim}.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{stratify.rtmle}},
#'   \code{\link{autoplot.stratified_rtmle}}
#' @examples
#' data(rtmle_object)
#' by_age <- stratify.rtmle(rtmle_object, age, breaks = c(-Inf, 50, Inf),
#'                         time_horizon = 1:2)
#' plot(by_age)
#' @method plot stratified_rtmle
#' @export
plot.stratified_rtmle <- function(x, ...) {
    autoplot.stratified_rtmle(x, ...)
}
