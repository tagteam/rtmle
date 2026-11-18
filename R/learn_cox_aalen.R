#' Learn censoring with a Cox-Aalen hazards model
#'
#' Fits \code{timereg::cox.aalen} using right-censored or start-stop data. A
#' formula such as \code{~ L + prop(A)} gives \code{L} an additive effect with
#' a time-varying coefficient and \code{A} a proportional effect. Covariates
#' are held fixed over each row's interval. With no \code{prop()} terms, the
#' model reduces to Aalen and is fitted with \code{timereg::aalen}.
#' Negative additive hazard jumps are truncated at zero with a warning.
#'
#' The learner currently rejects tied censoring times across subjects, because
#' \code{timereg} otherwise randomly breaks these ties. Cox censoring models
#' support tied censoring times.
#' Formulas with random effects are rejected. Additive predictors must be
#' identifiable in every censoring risk set; otherwise simplify the formula.
#' @inheritParams learn_aalen
#' @param ... Additional arguments for \code{timereg::cox.aalen}; by default,
#'   no simulation-based tests or robust variance calculations are requested.
#' @return A learner result with \code{predicted_values}, \code{fit}, and
#'   \code{fit_summary}. It can be supplied as \code{reuse_fit} for prediction.
#' @export
learn_cox_aalen <- function(character_formula, data, intervened_data = data,
                            time_variable = "time", time_ref = NULL,
                            save_fitted_objects = FALSE, reuse_fit = NULL,
                            id_variable = "id", left_limit = TRUE, ...) {
    result <- if (is.null(reuse_fit))
        continuous_timereg_fit(character_formula, data, "cox.aalen", id_variable, list(...)) else reuse_fit
    result$predicted_values <- continuous_timereg_survival(result, intervened_data,
                                                           time_variable, time_ref, left_limit)
    result$fit_summary <- if (is.null(result$fit)) numeric() else list(
        cumulative_coefficients = result$fit$cum, coefficients = result$fit$gamma)
    result
}
