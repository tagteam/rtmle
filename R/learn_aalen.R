#' Learn censoring with an Aalen additive hazards model
#'
#' Fits \code{timereg::aalen} using right-censored or start-stop data. The
#' covariates on each row are held fixed over its interval. Ordinary terms have
#' time-varying coefficients; \code{const()} specifies constant additive effects.
#' Negative fitted hazard jumps and negative constant hazard components are
#' truncated at zero with a warning when predicting survival.
#'
#' The learner currently rejects tied censoring times across subjects, because
#' \code{timereg} otherwise randomly breaks these ties. Cox censoring models
#' support tied censoring times.
#' Formulas with random effects are rejected. Additive predictors must be
#' identifiable in every censoring risk set; otherwise simplify the formula.
#' @inheritParams learn_coxph
#' @param character_formula A survival regression formula, supplied as a string
#'   or formula object.
#' @param id_variable Subject identifier for start-stop data. Defaults to \code{"id"}.
#' @param left_limit Whether to exclude a hazard jump at the interval endpoint.
#' @param ... Additional arguments for \code{timereg::aalen}; by default, no
#'   simulation-based tests or robust variance calculations are requested.
#' @details Negative fitted additive jumps and negative \code{const()} rates
#'   are truncated at zero with a warning so censoring survival stays in [0, 1].
#'   Nonzero offsets are not supported.
#' @return A learner result with \code{predicted_values}, \code{fit}, and
#'   \code{fit_summary}. It can be supplied as \code{reuse_fit} for prediction.
#' @export
learn_aalen <- function(character_formula, data, intervened_data = data,
                        time_variable = "time", time_ref = NULL,
                        save_fitted_objects = FALSE, reuse_fit = NULL,
                        id_variable = "id", left_limit = TRUE, ...) {
    result <- if (is.null(reuse_fit))
        continuous_timereg_fit(character_formula, data, "aalen", id_variable, list(...)) else reuse_fit
    result$predicted_values <- continuous_timereg_survival(result, intervened_data,
                                                           time_variable, time_ref, left_limit)
    result$fit_summary <- if (is.null(result$fit)) numeric() else list(
        cumulative_coefficients = result$fit$cum, coefficients = result$fit$gamma)
    result
}
