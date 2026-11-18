#' Learn continuous-time censoring with Weibull survival regression
#'
#' Fits a Weibull accelerated-failure-time model using
#' \code{survival::survreg}. Predictions are conditional survival probabilities
#' from \code{time_ref} through \code{time_variable}. The Weibull distribution
#' is continuous, so the left and right limits coincide. A model without any
#' observed censoring events predicts one. Stratified scale parameters are not
#' supported.
#' @param character_formula A survival regression formula as a string.
#' @param data Training data.
#' @param intervened_data Prediction data; defaults to \code{data}.
#' @param time_variable Column containing interval endpoints.
#' @param time_ref Optional column containing interval starting times.
#' @param save_fitted_objects Retained for learner-interface compatibility.
#' @param reuse_fit Result from a previous call to reuse for prediction.
#' @param ... Additional arguments passed to \code{survival::survreg}.
#' @return A list containing \code{predicted_values}, \code{fit}, and
#'   \code{fit_summary}.
#' @seealso \code{\link{learn_coxph}}, \code{\link{run_rtmle}}
#' @export
learn_weibull <- function(character_formula, data, intervened_data = data,
                          time_variable = "time", time_ref = NULL,
                          save_fitted_objects = FALSE, reuse_fit = NULL, ...) {
    if (is.null(reuse_fit)) {
        f <- stats::as.formula(character_formula)
        response <- stats::model.response(stats::model.frame(f, data, na.action = stats::na.fail))
        if (!inherits(response, "Surv") || ncol(response) != 2L ||
            !identical(attr(response, "type"), "right"))
            stop("learn_weibull requires a right-censored Surv(time, status) response.")
        if (any(response[, 1L] <= 0)) stop("Weibull training times must be positive.")
        if (!any(response[, 2L] == 1)) {
            fit <- NULL
        } else {
            fit <- survival::survreg(f, data = data, dist = "weibull", ...)
            if (length(fit$scale) != 1L) stop("Weibull stratified scale parameters are not supported.")
            fit$coefficients[is.na(fit$coefficients)] <- 0
        }
        result <- list(fit = fit)
    } else result <- reuse_fit

    endpoint <- intervened_data[[time_variable]]
    start <- if (is.null(time_ref)) rep(0, nrow(intervened_data)) else intervened_data[[time_ref]]
    if (!is.numeric(endpoint) || length(endpoint) != nrow(intervened_data) ||
        !is.numeric(start) || length(start) != nrow(intervened_data) ||
        any(!is.finite(endpoint)) || any(!is.finite(start)) || any(start < 0 | endpoint < start))
        stop("Censoring interval start and endpoint must be finite, nonnegative, and ordered.")
    if (is.null(result$fit)) {
        probabilities <- rep(1, nrow(intervened_data))
    } else {
        location <- stats::predict(result$fit, newdata = intervened_data, type = "lp")
        shape <- 1 / result$fit$scale
        # survreg uses a log-time location and scale. The ratio S(end)/S(start)
        # gives the probability of remaining uncensored over this interval.
        cumulative <- function(time) exp(shape * (log(time) - location))
        probabilities <- exp(-pmax(0, cumulative(endpoint) - cumulative(start)))
    }
    result$predicted_values <- continuous_nuisance_predictions(
        list(predicted_values = probabilities), nrow(intervened_data), "Censoring", TRUE)
    result$fit_summary <- if (is.null(result$fit)) numeric() else stats::coef(result$fit)
    result
}
