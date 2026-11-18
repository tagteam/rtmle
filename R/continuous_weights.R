continuous_weights <- function(data, a, last_event, regime = NULL) {
    n <- nrow(data)
    weights <- matrix(0, nrow = n, ncol = last_event)
    current <- rep(1, n)
    # At baseline and each subsequent event, carry forward the inverse
    # probability weight only for subjects still following the regime.
    for (j in 0:(last_event - 1)) {
        probability <- data[[paste0("propensity_", j)]]
        probability[is.na(probability)] <- 1 # No treatment decision at L events.
        if (j > 0) probability <- probability * data[[paste0("survival_censoring_", j)]]
        if (j == 0) {
            adheres <- rep(TRUE, n)
            decision_rows <- seq_len(n)
        } else {
            event <- data[[paste0("event_", j)]]
            adheres <- event %in% c("A", "L")
            decision_rows <- which(event == "A")
        }
        desired <- if (is.null(regime)) {
            rep(a, length(decision_rows))
        } else {
            continuous_intervention_values(data[decision_rows], regime, j)
        }
        adheres[decision_rows] <- data[[paste0("A_", j)]][decision_rows] == desired
        active <- which(current > 0 & !is.na(adheres) & adheres)
        updated <- numeric(n)
        updated[active] <- current[active] / probability[active]
        if (any(!is.finite(updated))) stop("Nonfinite weights: zero treatment or censoring probability on an adherent history.")
        weights[, j + 1] <- current <- updated
    }
    if (!any(weights[, 1] > 0)) stop("No baseline treatment support for the requested regime.")
    weights
}
