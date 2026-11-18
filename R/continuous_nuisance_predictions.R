continuous_nuisance_predictions <- function(result, n, what, positive = FALSE) {
    predictions <- as.numeric(result$predicted_values)
    if (length(predictions) != n || any(!is.finite(predictions)) ||
        any(predictions < 0 | predictions > 1) ||
        (positive && any(predictions == 0))) {
        stop(what, " predictions must contain one finite ",
             if (positive) "strictly positive " else "", "probability per prediction row.")
    }
    predictions
}
