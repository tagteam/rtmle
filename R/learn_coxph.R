#' Learn continuous-time censoring with a Cox model
#'
#' Uses the same learner arguments as \code{\link{learn_glm}}. Predictions are
#' conditional survival probabilities over an interval, using the left limit of
#' the cumulative hazard at its endpoint and the right limit at its start.
#' The formula may include \code{strata()} for stratum-specific hazards.
#' A \code{Surv(start, stop, status)} response fits a single model to
#' covariates that are constant on each recorded interval.
#' Frailty and random-effects formulas are not supported.
#'
#' @param character_formula A Cox regression formula, supplied as a string.
#' @param data Learning data.
#' @param intervened_data Prediction data. Defaults to \code{data}.
#' @param time_variable Column containing interval endpoints.
#' @param time_ref Optional column containing interval starting times.
#' @param save_fitted_objects Logical; retained for compatibility with other learners.
#' @param reuse_fit Optional result from a previous call.
#' @param left_limit Whether to exclude a cumulative-hazard jump at the
#'   prediction endpoint. Use \code{FALSE} at an interior interval boundary
#'   when combining consecutive intervals.
#' @param ... Arguments passed to \code{survival::coxph}.
#' @return A list containing \code{predicted_values}, the Cox \code{fit}, and
#'   \code{baseline_hazard}. With no censoring events, predictions equal one.
#' @export
learn_coxph <- function(character_formula, data, intervened_data = data,
                        time_variable = "time", time_ref = NULL,
                        save_fitted_objects = FALSE, reuse_fit = NULL,
                        left_limit = TRUE, ...) {
    if (!is.logical(left_limit) || length(left_limit) != 1L || is.na(left_limit))
        stop("left_limit must be TRUE or FALSE.")
    if (is.null(reuse_fit)) {
        f <- stats::as.formula(character_formula)
        # A covariate named "frailty" is still an ordinary predictor.
        calls <- sub("^.*::", "", continuous_formula_calls(f))
        if (any(grepl("^frailty($|\\.)", calls)) || any(calls %in% c("|", "||")))
            stop("Censoring models do not support frailty or random effects.")
        # Make the usual unqualified strata() usable without attaching survival.
        formula_env <- new.env(parent = environment(f))
        formula_env$strata <- survival::strata
        environment(f) <- formula_env
        response <- stats::model.response(stats::model.frame(f, data, na.action = stats::na.fail))
        if (!inherits(response, "Surv") ||
            !attr(response, "type") %in% c("right", "counting")) {
            stop("learn_coxph requires a right-censored Surv(time, status) or Surv(start, stop, status) response.")
        }
        if (any(!is.finite(response)) || any(response[, 1L] < 0) ||
            (ncol(response) == 3L && any(response[, 2L] <= response[, 1L])))
            stop("Cox training times must be finite, nonnegative, and ordered.")
        if (!any(response[, ncol(response)] == 1)) {
            result <- list(fit = NULL, baseline_hazard = data.frame(time = numeric(), hazard = numeric()))
        } else {
            fit <- survival::coxph(f, data = data, x = TRUE, ...)
            # Aliased coefficients have zero contribution to predictions.
            fit$coefficients[is.na(fit$coefficients)] <- 0
            result <- list(fit = fit, baseline_hazard = survival::basehaz(fit, centered = FALSE))
        }
    } else {
        result <- reuse_fit
    }
    result$predicted_values <- continuous_cox_survival(result, intervened_data,
                                                      time_variable, time_ref, left_limit)
    result$fit_summary <- if (!is.null(result$fit)) stats::coef(result$fit) else numeric()
    result
}
