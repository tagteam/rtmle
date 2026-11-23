#' Print stratified rtmle analyses
#'
#' Prints the number of participants and the estimates for each stratum.
#'
#' @param x A fitted \code{stratified_rtmle} object.
#' @param ... Arguments passed to \code{\link{summary.stratified_rtmle}}, such
#'   as \code{targets}, \code{reference}, and \code{digits}.
#' @return The input object, invisibly.
#' @seealso \code{\link{stratify.rtmle}},
#'   \code{\link{summary.stratified_rtmle}}
#' @examples
#' data(rtmle_object)
#' by_sex <- stratify.rtmle(rtmle_object, sex, time_horizon = 1:2)
#' print(by_sex)
#' @method print stratified_rtmle
#' @export
print.stratified_rtmle <- function(x, ...) {
    variable <- attr(x, "variable", exact = TRUE)
    cat("RTMLE analyses stratified by ", variable, "\n\n", sep = "")
    counts <- data.frame(Stratum = names(x),
                         Participants = unname(attr(x, "counts", exact = TRUE)))
    count_column <- if (variable == "Participants") "Participant count" else "Participants"
    names(counts) <- c(variable, count_column)
    print(counts, row.names = FALSE)
    cat("\nEstimates\n")
    estimates <- summary(x, ...)
    columns <- c(names(estimates)[1L], "Target", "Regime", "Reference", "Target_parameter",
                 "Time_horizon", "Estimator", "Estimate (CI_95)")
    print(estimates[, intersect(columns, names(estimates)), with = FALSE])
    invisible(x)
}

