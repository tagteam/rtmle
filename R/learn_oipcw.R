#' Learner for sequential regressions for continuous-time analyses
#'
#' Fits the logistic estimating equation \eqn{X^T(Y - expit(X\beta))=0}.
#' The non-negative response may exceed one, as happens with inverse censoring
#' weighting. Predictions remain in \eqn{[0,1]}. An unsolved estimating equation
#' raises an error rather than returning an unverified fit.
#' The method is strongly inspired by Scheike, Zhang, and Gerds (2008),
#' i.e., https://academic.oup.com/biomet/article-abstract/95/1/205/219426
#'
#' @param character_formula Regression formula as a character string.
#' @param data Data used for learning.
#' @param intervened_data Data used for prediction.
#' @param save_fitted_objects Whether to return the fitted model as \code{fit}.
#' @param reuse_fit A previous learner result containing a saved \code{fit}.
#' @param maxit Maximum number of Newton iterations.
#' @param tol Tolerance for the estimating equation.
#' @param ... Additional arguments, currently unused.
#' @return A list containing \code{predicted_values}, \code{fit_summary}, and,
#'   when requested, \code{fit}.
#' @seealso \code{\link{run_rtmle}}, \code{\link{learn_glm}}
#' @export
learn_oipcw <- function(character_formula, data, intervened_data,
                        save_fitted_objects = FALSE, reuse_fit = NULL,
                        maxit = 1000, tol = 1e-8, ...) {
    if (is.null(reuse_fit)) {
        frame <- stats::model.frame(stats::formula(character_formula), data,
                                     na.action = stats::na.fail)
        Y <- stats::model.response(frame)
        if (!is.numeric(Y) || !length(Y) || any(!is.finite(Y)) || any(Y < 0)) {
            stop("learn_oipcw requires finite, non-negative numeric responses.")
        }
        Y <- as.numeric(Y)
        # A constant response in [0, 1] has a constant fit, including the
        # boundary cases where a finite logistic intercept does not exist.
        if (all(Y == Y[1]) && Y[1] <= 1) {
            fit <- list(constant = Y[1],
                        coefficients = c("(Intercept)" = stats::qlogis(Y[1])))
        } else {
            terms <- stats::terms(frame)
            # Single-level factors cannot be encoded by model.matrix. Remove
            # their terms; numeric constant and aliased columns are handled by QR.
            constant_factors <- names(frame)[vapply(frame, function(x) {
                is.factor(x) && length(unique(x)) == 1L
            }, logical(1))]
            if (length(constant_factors)) {
                factors <- attr(terms, "factors")
                removed <- which(colSums(factors[constant_factors, , drop = FALSE]) > 0)
                terms <- stats::drop.terms(terms, removed, keep.response = TRUE)
                frame <- stats::model.frame(terms, data, na.action = stats::na.fail)
            }
            X <- stats::model.matrix(terms, frame)
            # Solve only for independent columns; restore zero coefficients
            # for aliased columns so prediction uses the original design.
            decomposition <- qr(X)
            if (!decomposition$rank) stop("The formula has no estimable coefficients.")
            columns <- decomposition$pivot[seq_len(decomposition$rank)]
            X_fit <- X[, columns, drop = FALSE]
            intercept <- match("(Intercept)", colnames(X_fit))
            mean_Y <- mean(Y)
            if (!is.na(intercept) && mean_Y >= 1) {
                stop("The IPCW logistic estimating equation has no finite solution: the outcome mean is at least one.")
            }
            beta <- numeric(ncol(X_fit))
            if (!is.na(intercept)) beta[intercept] <- stats::qlogis(mean_Y)
            offset <- stats::model.offset(frame)
            if (is.null(offset)) offset <- rep(0, length(Y))
            # Unlike a binomial likelihood, this estimating equation allows
            # IPCW pseudo-outcomes above one while keeping predictions bounded.
            beta <- as.vector(estimating_equation_cpp(
                X = X_fit, Y = Y, beta = beta, offset = offset,
                maxit = maxit, tol = tol
            ))
            coefficients <- stats::setNames(numeric(ncol(X)), colnames(X))
            coefficients[columns] <- beta
            fit <- list(terms = stats::delete.response(terms),
                        coefficients = coefficients,
                        contrasts = attr(X, "contrasts"),
                        xlevels = lapply(frame[vapply(frame, is.factor, logical(1))], levels))
        }
    } else {
        fit <- reuse_fit$fit
        if (is.null(fit)) stop("reuse_fit must contain a saved fit.")
    }
    if (!is.null(fit$constant)) {
        predicted_values <- rep(fit$constant, NROW(intervened_data))
    } else {
        prediction_frame <- stats::model.frame(fit$terms, intervened_data,
                                                xlev = fit$xlevels,
                                                na.action = stats::na.fail)
        X_new <- stats::model.matrix(fit$terms, prediction_frame,
                                      contrasts.arg = fit$contrasts)
        offset <- stats::model.offset(prediction_frame)
        if (is.null(offset)) offset <- rep(0, nrow(X_new))
        predicted_values <- as.vector(stats::plogis(
            X_new[, names(fit$coefficients), drop = FALSE] %*% fit$coefficients + offset
        ))
    }
    output <- list(predicted_values = predicted_values,
                   fit_summary = matrix(fit$coefficients, ncol = 1,
                                        dimnames = list(names(fit$coefficients), "Estimate")))
    if (isTRUE(save_fitted_objects)) output$fit <- fit
    output
}
