continuous_sequential_regression <- function(nuisance, th, a, learner, estimator,
                                              lag, seed, prediction_range,
                                              scale_pseudo_outcome, stratify,
                                              intervene_all, intervention_hook,
                                              save_fitted_objects, regime = NULL,
                                              model_options = NULL) {
    prepared <- nuisance$prepared_data_object
    last_event <- prepared$info$last_event[match(th, prepared$info$time_horizon)]
    data <- data.table::copy(nuisance$data)
    # The final pooled event is specific to this horizon. Earlier events and
    # their nuisance predictions are shared across the horizon fits.
    for (variable in c("time", "event", "survival_censoring")) {
        column <- paste0(variable, "_", last_event)
        data.table::set(data, j = column, value = data[[paste0(variable, "_pooled_", last_event)]])
    }
    weights <- continuous_weights(data, a, last_event, regime)
    n <- nrow(data)
    next_Q <- residuals <- ipw <- numeric(n)
    models <- diagnostics <- list()
    dynamic <- isTRUE(regime$dynamic_intervention)
    outcome_name <- model_options$node_names$outcome
    if (is.null(outcome_name) || !length(outcome_name)) outcome_name <- "Y"
    # Work backwards from the horizon to baseline. next_Q contains the
    # intervened prediction from the next event for each subject.
    for (k in rev(seq_len(last_event))) {
        previous_event <- data[[paste0("event_", k - 1)]]
        previous_time <- data[[paste0("time_", k - 1)]]
        rows <- which(previous_event %in% c("A", "L") & previous_time < th)
        if (!length(rows)) next
        training <- data.table::copy(data[rows])
        event <- training[[paste0("event_", k)]]
        before_horizon <- training[[paste0("time_", k)]] <= th
        # The outcome is one at an observed outcome event, the next Q at a
        # treatment or covariate event, and zero after other terminal events
        # or beyond the horizon. Divide by survival over the interval for IPCW.
        response <- as.numeric(before_horizon & event == "Y") +
            as.numeric(before_horizon & event %in% c("A", "L")) * next_Q[rows]
        G <- training[[paste0("survival_censoring_", k)]]
        if (any(!is.finite(G)) || any(G <= 0 | G > 1)) stop("Invalid censoring survival probabilities.")
        response <- response / G
        data.table::set(training, j = "rtmle_pseudo_outcome", value = response)
        # Fit on observed histories and predict on histories with treatment set
        # by the regime. By default only the most recent decision is changed.
        prediction <- data.table::copy(training)
        nodes <- if (intervene_all) 0:(k - 1) else k - 1
        for (j in nodes) {
            if (j > 0 && intervene_all && dynamic) {
                measurement_rows <- which(prediction[[paste0("event_", j)]] == "L")
                data.table::set(prediction, i = measurement_rows, j = paste0("A_", j),
                                value = prediction[[paste0("A_", j - 1)]][measurement_rows])
            }
            intervention_rows <- if (j == 0 || (intervene_all && !dynamic)) {
                seq_len(nrow(prediction))
            } else {
                which(prediction[[paste0("event_", j)]] == "A")
            }
            desired <- if (is.null(regime)) {
                rep(a, length(intervention_rows))
            } else {
                continuous_intervention_values(prediction[intervention_rows], regime, j)
            }
            data.table::set(prediction, i = intervention_rows, j = paste0("A_", j), value = desired)
        }
        if (!is.null(intervention_hook)) {
            prediction <- intervention_hook(prediction, k - 1)
            if (!is.data.frame(prediction) || !identical(prediction$id, training$id))
                stop("intervention_hook must preserve prediction rows and subject identifiers.")
            prediction <- data.table::as.data.table(prediction)
        }
        learning_rows <- if (stratify) which(weights[rows, k] > 0) else seq_along(rows)
        if (!length(learning_rows)) stop("No adherent subjects available for stratified regression at event ", k, ".")
        covariates <- continuous_history(training[learning_rows], prepared$time_covariates,
                                          prepared$baseline_covariates, "pseudo_outcome", lag, k,
                                          model_options = model_options,
                                          outcome_name = paste0(outcome_name, "_", k))
        formula <- continuous_nuisance_formula("rtmle_pseudo_outcome", covariates)
        # Ordinary bounded outcome learners need a response in [0, 1]. Undo
        # this scaling immediately after prediction; the IPCW outcome itself
        # can exceed one.
        response_scale <- if (scale_pseudo_outcome) max(1, response[learning_rows]) else 1
        data.table::set(training, j = "rtmle_pseudo_outcome", value = response / response_scale)
        fitted <- fitter(intervention_node = k, learner = learner, formula = formula,
                         data = training[learning_rows], intervened_data = prediction,
                         id_variable = "id", minority_threshold = 0, seed = seed,
                         save_fitted_objects = save_fitted_objects)
        diagnostics <- continuous_merge_diagnostics(diagnostics, fitted$diagnostics)
        Q <- as.numeric(fitted$predicted_values) * response_scale
        if (length(Q) != length(rows) || any(!is.finite(Q))) stop("Outcome learner must return one finite prediction per row.")
        Q <- pmax(0, pmin(1, Q))
        epsilon <- NULL
        if (estimator == "tmle") {
            update <- continuous_tmle_update(response, Q, weights[rows, k], prediction_range)
            Q <- update$predicted_values
            epsilon <- update$epsilon
        }
        # Accumulate the influence-curve residual and the direct IPW estimate
        # before carrying the intervened Q backwards to the preceding event.
        residuals[rows] <- residuals[rows] + weights[rows, k] * (response - Q)
        ipw[rows] <- ipw[rows] + weights[rows, k] * as.numeric(before_horizon & event == "Y") / G
        next_Q[rows] <- Q
        # Retain summaries; full objects follow run_rtmle's existing setting.
        models[[paste0("event_", k)]] <- list(formula = formula,
                                              fit_summary = fitted$fit_summary,
                                              fit = if (save_fitted_objects) fitted$fit,
                                              response_scale = response_scale,
                                              epsilon = epsilon)
    }
    # After the baseline iteration next_Q is each subject's baseline Q.
    ice_ipcw <- mean(next_Q)
    ic <- residuals + next_Q - ice_ipcw
    estimate <- if (estimator == "one-step") ice_ipcw + mean(ic) else ice_ipcw
    # A plug-in g-formula analysis has no influence-based standard error here.
    se <- if (estimator == "g-formula") NA_real_ else stats::sd(ic) / sqrt(n)
    list(estimate = estimate, se = se, ic = if (estimator == "g-formula") rep(NA_real_, n) else ic,
         ice_ipcw = if (estimator == "tmle") NA_real_ else ice_ipcw,
         ipw = mean(ipw), models = models, diagnostics = diagnostics)
}
