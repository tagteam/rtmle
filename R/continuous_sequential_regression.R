continuous_sequential_regression <- function(nuisance, th, a, learner, estimator,
                                              lag, seed, prediction_range,
                                              stratify,
                                              intervene_all, intervention_hook,
                                              save_fitted_objects, regime = NULL,
                                              model_options = NULL, q_grid = NULL,
                                             grid_learner = learner) {
    prepared <- nuisance$prepared_data_object
    last_event <- prepared$info$last_event[match(th, prepared$info$time_horizon)]
    data <- data.table::copy(nuisance$data)
    # The final pooled event is specific to this horizon. Earlier events and
    # their nuisance predictions are shared across the horizon fits.
    for (variable in c("time", "event", "survival_censoring")) {
        column <- paste0(variable, "_", last_event)
        data.table::set(data, j = column, value = data[[paste0(variable, "_pooled_", last_event)]])
    }
    fit_propensity <- estimator != "sequential_regression"
    weights <- if (fit_propensity || stratify)
        continuous_weights(data, a, last_event, regime, inverse_probability = fit_propensity)
    n <- nrow(data)
    next_Q <- numeric(n)
    if (fit_propensity) residuals <- ipw <- numeric(n)
    models <- diagnostics <- list()
    if (!is.null(q_grid)) {
        if (!is.numeric(q_grid) || length(q_grid) < 2L || any(!is.finite(q_grid)) ||
            q_grid[1L] != 0 || q_grid[length(q_grid)] != th || any(diff(q_grid) <= 0)) {
            stop("q_grid must increase from zero to the analysis horizon.")
        }
        grid_predictions <- vector("list", last_event)
    }
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
        event_time <- training[[paste0("time_", k)]]
        before_horizon <- event_time <= th
        # The outcome is one at an observed outcome event, the next Q at a
        # treatment or covariate event, and zero after other terminal events
        # or beyond the horizon. Divide by survival over the interval for IPCW.
        unweighted <- as.numeric(before_horizon & event == "Y") +
            as.numeric(before_horizon & event %in% c("A", "L")) * next_Q[rows]
        G <- training[[paste0("survival_censoring_", k)]]
        if (any(!is.finite(G)) || any(G <= 0 | G > 1)) stop("Invalid censoring survival probabilities.")
        response <- unweighted / G
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
        if (!fit_propensity && !stratify && k == 1L &&
            !any(training$A_0 == prediction$A_0, na.rm = TRUE)) {
            stop("No baseline treatment support for the requested regime.")
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
        fit_q <- function(cutoff, train_rows, prediction_data, save_fit, outcome_learner) {
            observed_by_cutoff <- event_time <= cutoff
            if (is_ipcw_glm_learner(outcome_learner)) {
                # Before the cutoff, only observed non-dropout events contribute.
                # Beyond it the pseudo-outcome is zero and observed at the cutoff.
                at_cutoff <- observed_by_cutoff[train_rows]
                survival <- G[train_rows]
                later <- which(!at_cutoff)
                if (length(later)) {
                    survival[later] <- continuous_censoring_at_cutoff(
                        nuisance, k, k == last_event, rows[train_rows[later]], cutoff)
                }
                ipcw <- as.numeric(!at_cutoff | event[train_rows] != "C") / survival
                fit_data <- training[train_rows]
                data.table::set(fit_data, j = "rtmle_pseudo_outcome", value = unweighted[train_rows] * at_cutoff)
                fit_covariates <- continuous_varying_covariates(fit_data, covariates)
                weighted_formula <- continuous_nuisance_formula("rtmle_pseudo_outcome", fit_covariates)
                warnings <- character()
                fit <- withCallingHandlers(do.call(outcome_learner$fun, c(list(
                    character_formula = weighted_formula, data = fit_data,
                    intervened_data = prediction_data, weights = ipcw,
                    save_fitted_objects = save_fit), outcome_learner$args)), warning = function(w) {
                    warnings <<- c(warnings, conditionMessage(w))
                    invokeRestart("muffleWarning")
                })
                if (length(warnings)) {
                    fit$diagnostics <- list(learner_warnings = data.table::data.table(
                        learner = outcome_learner$name, warning = warnings))
                }
            } else {
                data.table::set(training, j = "rtmle_pseudo_outcome", value = response * observed_by_cutoff)
                fit <- fitter(intervention_node = k, learner = outcome_learner, formula = formula,
                              data = training[train_rows], intervened_data = prediction_data,
                              id_variable = "id", minority_threshold = 0, seed = seed,
                              save_fitted_objects = save_fit)
            }
            values <- as.numeric(fit$predicted_values)
            if (length(values) != nrow(prediction_data) || any(!is.finite(values)))
                stop("Outcome learner must return one finite prediction per row.")
            fit$predicted_values <- pmax(0, pmin(1, values))
            fit
        }
        fitted <- fit_q(th, learning_rows, prediction, save_fitted_objects, learner)
        diagnostics <- continuous_merge_diagnostics(diagnostics, fitted$diagnostics)
        Q <- fitted$predicted_values
        epsilon <- NULL
        if (estimator == "tmle") {
            update <- continuous_tmle_update(response, Q, weights[rows, k], prediction_range)
            Q <- update$predicted_values
            epsilon <- update$epsilon
        }
        if (!is.null(q_grid)) {
            # q(u) truncates the time of this interval's event. The continuation
            # value at a treatment or covariate event is still evaluated at th.
            values <- matrix(0, nrow = n, ncol = length(q_grid))
            values[rows, length(q_grid)] <- Q
            entry <- previous_time[rows]
            for (j in which(q_grid > 0 & q_grid < th)) {
                u <- q_grid[j]
                eligible <- which(entry < u)
                train_rows <- learning_rows[entry[learning_rows] < u]
                if (!any(event_time[train_rows] <= u & unweighted[train_rows] > 0)) next
                # Grid refinements must not change the random draws used by
                # stochastic learners for the original horizon regressions.
                grid_fit <- tryCatch(
                    withr::with_preserve_seed(fit_q(u, train_rows, prediction[eligible], FALSE,
                                                    grid_learner)),
                    error = function(e) stop("Grid outcome regression failed at event ", k,
                        ", u = ", format(u), " (tau = ", format(th), "): ",
                        conditionMessage(e), call. = FALSE)
                )
                values[rows[eligible], j] <- grid_fit$predicted_values
                diagnostics <- continuous_merge_diagnostics(diagnostics, grid_fit$diagnostics)
            }
            grid_predictions[[k]] <- values
        }
        # Accumulate the influence-curve residual and the direct IPW estimate
        # before carrying the intervened Q backwards to the preceding event.
        if (fit_propensity) {
            residuals[rows] <- residuals[rows] + weights[rows, k] * (response - Q)
            ipw[rows] <- ipw[rows] + weights[rows, k] * as.numeric(before_horizon & event == "Y") / G
        }
        next_Q[rows] <- Q
        # Retain summaries; full objects follow run_rtmle's existing setting.
        models[[paste0("event_", k)]] <- list(formula = formula,
                                              fit_summary = fitted$fit_summary,
                                              fit = if (save_fitted_objects) fitted$fit,
                                              epsilon = epsilon)
    }
    # After the baseline iteration next_Q is each subject's baseline Q.
    ice_ipcw <- mean(next_Q)
    ic <- if (fit_propensity) residuals + next_Q - ice_ipcw else rep(NA_real_, n)
    estimate <- if (estimator == "one-step") ice_ipcw + mean(ic) else ice_ipcw
    se <- if (fit_propensity) stats::sd(ic) / sqrt(n) else NA_real_
    result <- list(estimate = estimate, se = se, ic = ic,
                   ice_ipcw = if (estimator == "tmle") NA_real_ else ice_ipcw,
                   ipw = if (fit_propensity) mean(ipw) else NA_real_,
                   models = models, diagnostics = diagnostics)
    if (!is.null(q_grid)) {
        result$q_grid <- grid_predictions
        result$weights <- weights
    }
    result
}
