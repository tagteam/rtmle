# Fit and apply one intercept fluctuation to continuous-time IPCW outcomes.
# Unlike a binomial GLM, this score also permits responses exceeding one.
continuous_tmle_update <- function(Y, Q, weights,
                                   prediction_range = c(0.0001, 0.9999)) {
    if (!length(Y) || length(Q) != length(Y) || length(weights) != length(Y)) {
        stop("Y, Q, and weights must have the same length and be non-empty.")
    }
    if (!is.numeric(Y) || !is.numeric(Q) || !is.numeric(weights) ||
        any(!is.finite(c(Y, Q, weights)))) {
        stop("Y, Q, and weights must contain only finite numeric values.")
    }
    if (any(Y < 0) || any(weights < 0)) {
        stop("Y and weights must be non-negative.")
    }
    if (any(Q < 0 | Q > 1)) stop("Q must lie in [0, 1].")
    if (!is.numeric(prediction_range) || length(prediction_range) != 2L ||
        any(!is.finite(prediction_range)) || prediction_range[1] <= 0 ||
        prediction_range[2] >= 1 || prediction_range[1] >= prediction_range[2]) {
        stop("prediction_range must contain two increasing values inside (0, 1).")
    }
    positive <- weights > 0
    if (!any(positive)) {
        return(list(predicted_values = Q, epsilon = 0, score = 0))
    }

    # Solve sum_i w_i (Y_i - expit(logit(Q_i) + epsilon)) = 0.
    # Normalizing positive weights does not change the root and avoids large
    # intermediate values. A root exists only if the weighted mean Y is in
    # [0, 1], because the updated predictions are bounded probabilities.
    normalized_weights <- weights[positive] / max(weights[positive])
    weight_sum <- sum(normalized_weights)
    target_mean <- sum(normalized_weights * Y[positive]) / weight_sum
    normalized_weights <- normalized_weights / weight_sum
    # Allow only floating-point roundoff at the upper boundary. In particular,
    # an all-one response must not fail because normalized weights sum above one.
    boundary_tolerance <- 10 * .Machine$double.eps
    if (target_mean > 1 + boundary_tolerance) {
        stop("The continuous TMLE score has no solution: the weighted outcome mean exceeds one.")
    }
    if (target_mean >= 1 - boundary_tolerance) target_mean <- 1
    if (target_mean == 0 || target_mean == 1) {
        # The score reaches zero only in the limit at these boundaries.
        predicted_values <- rep(target_mean, length(Q))
        epsilon <- if (target_mean == 0) -Inf else Inf
    } else {
        bounded_Q <- pmax(prediction_range[1], pmin(prediction_range[2], Q))
        offset <- stats::qlogis(bounded_Q)
        score <- function(epsilon) {
            sum(normalized_weights *
                    (Y[positive] - stats::plogis(offset[positive] + epsilon)))
        }
        epsilon <- stats::uniroot(score, interval = c(-50, 50),
                                  extendInt = "downX", tol = 1e-12)$root
        predicted_values <- stats::plogis(offset + epsilon)
    }
    list(predicted_values = predicted_values, epsilon = epsilon,
         score = sum(weights * (Y - predicted_values)))
}
