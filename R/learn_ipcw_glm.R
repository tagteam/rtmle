#' Learner for sequential regressions for continuous-time analyses
#'
#' Fits a quasi-binomial GLM to the unweighted pseudo-outcome in \eqn{[0,1]},
#' with inverse censoring probabilities supplied as observation weights.
#' The method is strongly inspired by Blanche et al. (2023), i.e.,
#' https://link.springer.com/article/10.1007/s10985-022-09564-6.
#'
#' @param character_formula Regression formula as a character string.
#' @param data Data used for learning.
#' @param intervened_data Data used for prediction.
#' @param weights Nonnegative observation weights, in the same order as \code{data}.
#' @param link Link for \code{stats::quasibinomial()}, defaulting to \code{"logit"}.
#' @param maxit Maximum number of GLM iterations.
#' @param save_fitted_objects Whether to return the fitted model as \code{fit}.
#' @param reuse_fit A previous learner result containing a saved \code{fit}.
#' @param ... Additional arguments, currently unused.
#' @return A list containing \code{predicted_values}, \code{fit_summary}, and,
#'   when requested, \code{fit}.
#' @seealso \code{\link{run_rtmle}}, \code{\link{learn_oipcw}}
#' @export
learn_ipcw_glm <- function(character_formula, data, intervened_data, weights,
                           link = "logit", maxit = 100, save_fitted_objects = FALSE,
                           reuse_fit = NULL, ...) {
    if (is.null(reuse_fit)) {
        if (!is.numeric(weights) || length(weights) != nrow(data) ||
            any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
            stop("learn_ipcw_glm requires one finite, nonnegative weight per row and at least one positive weight.")
        }
        formula <- stats::as.formula(character_formula)
        .rtmle_ipcw_weights <- weights
        fit <- stats::glm(formula, data = data, weights = .rtmle_ipcw_weights,
                          family = stats::quasibinomial(link = link),
                          control = stats::glm.control(maxit = maxit), na.action = stats::na.fail)
    } else {
        fit <- reuse_fit$fit
        if (is.null(fit)) stop("reuse_fit must contain a saved fit.")
    }
    coefficients <- stats::coef(fit)
    result <- list(predicted_values = as.numeric(stats::predict(fit, intervened_data, type = "response")),
                   fit_summary = matrix(coefficients, ncol = 1L,
                                        dimnames = list(names(coefficients), "Estimate")))
    if (isTRUE(save_fitted_objects)) result$fit <- fit
    result
}
