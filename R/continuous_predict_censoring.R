# Integrate along the recorded intervals. Pooled endpoints may cross several
# updates. Only the final endpoint uses a left limit; interior jumps count once.
continuous_predict_censoring <- function(model, prediction, endpoint, start, specification, call_learner) {
    requests <- data.table::data.table(id = prediction$id,
        .rtmle_from = prediction[[start]], .rtmle_to = prediction[[endpoint]],
        .rtmle_request = seq_len(nrow(prediction)))
    if (anyNA(requests) || any(!is.finite(requests$.rtmle_from)) ||
        any(!is.finite(requests$.rtmle_to)) || any(requests$.rtmle_to < requests$.rtmle_from))
        stop("Invalid censoring prediction intervals.")
    pieces <- merge(requests, model$intervals, by = "id", allow.cartesian = TRUE)
    keep <- pieces$.rtmle_stop > pieces$.rtmle_from & pieces$.rtmle_start < pieces$.rtmle_to
    pieces <- pieces[keep]
    probabilities <- rep(1, nrow(prediction))
    data.table::set(pieces, j = ".rtmle_start", value = pmax(pieces$.rtmle_start, pieces$.rtmle_from))
    data.table::set(pieces, j = ".rtmle_stop", value = pmin(pieces$.rtmle_stop, pieces$.rtmle_to))
    covered <- numeric(nrow(prediction))
    if (nrow(pieces)) {
        durations <- rowsum(pieces$.rtmle_stop - pieces$.rtmle_start, pieces$.rtmle_request)
        covered[as.integer(rownames(durations))] <- durations[, 1L]
    }
    if (any(abs(covered - (requests$.rtmle_to - requests$.rtmle_from)) >
            1e-8 * pmax(1, covered)))
        stop("Censoring prediction extends beyond the recorded event history.")
    if (!nrow(pieces)) return(probabilities)
    final <- pieces$.rtmle_stop == pieces$.rtmle_to
    values <- numeric(nrow(pieces))
    for (left in c(FALSE, TRUE)) {
        rows <- which(final == left)
        if (!length(rows)) next
        result <- call_learner(c(list(character_formula = model$formula, data = model$intervals,
            intervened_data = pieces[rows], time_variable = ".rtmle_stop", time_ref = ".rtmle_start",
            reuse_fit = model$fitted, left_limit = left, save_fitted_objects = TRUE),
            specification$learner$args))
        values[rows] <- continuous_nuisance_predictions(result, length(rows), "Censoring", TRUE)
    }
    totals <- rowsum(-log(values), pieces$.rtmle_request, reorder = FALSE)
    probabilities[as.integer(rownames(totals))] <- exp(-totals[, 1L])
    probabilities
}
