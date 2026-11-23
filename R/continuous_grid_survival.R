# S_k is the survival of the next full-data event, conditional on the previous
# history. Each event index has its own gap-time model; A, L, Y and D are
# failures, while dropout and the administrative end right-censor that time.
continuous_event_specification <- function(learner) {
    parsed <- continuous_nuisance_learner(learner)
    if (identical(parsed$name, "superlearn"))
        stop("event_learner must be a single survival learner.")
    for (option in c("marginal_censoring", "time_varying")) {
        setting <- parsed$args[[option, exact = TRUE]]
        if (!is.null(setting) && !identical(setting, FALSE))
            stop("event_learner requires event-specific models conditional on prior history; ",
                 option, " must be FALSE.")
        parsed$args[[option]] <- NULL
    }
    if (any(c("marginal_censoring_formula", "censoring_formula") %in% names(parsed$args)))
        stop("Censoring formula options do not apply to event_learner.")
    parsed
}

# Input is a private copy of the risk set. Express both the response time
# and the preceding event times relative to entry into this interval.
continuous_gap_history <- function(data, k, endpoint) {
    start <- data[[paste0("time_", k - 1L)]]
    data.table::set(data, j = endpoint, value = data[[endpoint]] - start)
    if (k > 1L) for (j in seq_len(k - 1L)) {
        column <- paste0("time_", j)
        data.table::set(data, j = column, value = start - data[[column]])
    }
    data
}

# Grid times u are on the absolute follow-up clock, and may differ by row.
# Use S_k(u-) at dropout jumps. For a quadrature cell (from, u], evaluate
# S_k just after 'from' for its compensator, and S^c_k just after u.
continuous_grid_survival <- function(nuisance, k, pooled, rows, u, from, jump_rows) {
    grid <- nuisance$grid_survival
    endpoint <- paste0(if (pooled) "time_pooled_" else "time_", k)
    start_name <- paste0("time_", k - 1L)
    start <- nuisance$data[[start_name]][rows]
    model <- grid$event[[endpoint]]
    index <- match(rows, model$rows)
    if (is.null(model) || anyNA(index) || length(u) != length(rows) ||
        any(!is.finite(u)) || any(u < start))
        stop("Invalid risk set or time for grid survival prediction.")
    predict_model <- function(fit, learner, training, newdata, time_ref, left, label, positive = TRUE) {
        result <- do.call(learner$fun, c(list(character_formula = fit$formula,
            data = training, intervened_data = newdata, time_variable = endpoint,
            time_ref = time_ref, left_limit = left, reuse_fit = fit$fit,
            save_fitted_objects = TRUE), learner$args))
        continuous_nuisance_predictions(result, nrow(newdata), label, positive)
    }
    predict_event <- function(at, positions, left) {
        # The conditional survival at interval entry is one. This also avoids
        # fitting a prediction matrix for an empty set of dropout jumps.
        S <- rep(1, length(positions))
        later <- which(at > start[positions])
        if (length(later)) {
            use <- positions[later]
            prediction <- model$training[index[use]]
            data.table::set(prediction, j = endpoint, value = at[later] - start[use])
            S[later] <- predict_model(model, grid$event_learner, model$training, prediction,
                                      NULL, left, "Event survival", positive = FALSE)
        }
        S
    }
    if (length(from) != length(rows) || any(!is.finite(from)) ||
        any(from < start | from > u))
        stop("Invalid cell start for grid survival prediction.")
    survival <- predict_event(from, seq_along(rows), FALSE)
    censoring <- grid$censoring
    if (!censoring$marginal && !censoring$time_varying) {
        conditional <- grid$conditional_censoring[[endpoint]]
        prediction <- conditional$training[index]
        data.table::set(prediction, j = endpoint, value = u - start)
        G <- predict_model(conditional, censoring$learner, conditional$training,
                           prediction, NULL, FALSE, "Censoring")
    } else {
        prediction <- nuisance$data[rows]
        data.table::set(prediction, j = endpoint, value = u)
        if (censoring$marginal) {
            G <- predict_model(grid$marginal_censoring, censoring$learner,
                               nuisance$prepared_data_object$data_marginal_censoring,
                               prediction, start_name, FALSE, "Censoring")
        } else {
            G <- continuous_predict_censoring(grid$shared_censoring, prediction, endpoint,
                start_name, censoring, function(args) do.call(censoring$learner$fun, args),
                left_limit = FALSE)
        }
    }
    list(event = survival, censoring = G,
         event_at_jump = predict_event(u[jump_rows], jump_rows, TRUE))
}
