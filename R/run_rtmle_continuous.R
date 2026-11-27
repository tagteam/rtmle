# Continuous-time ICE-IPCW recursion. The models and results use the same
# learner protocol and result tables as the discrete-time analysis.
run_rtmle_continuous <- function(x, targets = NULL, learner = NULL,
                                 estimator = "tmle", time_horizon = NULL,
                                 refit = TRUE, seed = NULL, subsets = NULL,
                                 keep_influence = TRUE, save_fitted_objects = FALSE,
                                 progressbar = 0, verbose = FALSE,
                                 propensity_learner = if (is.null(learner)) "learn_glm" else learner,
                                 censoring_learner = "learn_coxph",
                                 conservative = TRUE, lag = NULL, gbound = NULL,
                                 stratify_sequential_regression = FALSE,
                                 intervene_all_sequential_regression = FALSE,
                                 intervention_hook = NULL,
                                 exclude_latest_covariate = NULL,
                                 event_learner = "learn_coxph",
                                 grid_learner = learner,
                                 integration_grid_size = 100L, ...) {
    if (!is.logical(conservative) || length(conservative) != 1L || is.na(conservative))
        stop("conservative must be TRUE or FALSE.")
    if (!is.character(estimator) || length(estimator) != 1L ||
        !estimator %in% c("tmle", "one-step", "sequential_regression", "g-formula"))
        stop("Continuous estimator must be 'tmle', 'one-step', or 'sequential_regression'.")
    if (estimator == "g-formula") estimator <- "sequential_regression"
    if (!conservative) {
        if (estimator != "one-step")
            stop("The experimental non-conservative correction requires estimator = 'one-step'.")
        if (is.null(event_learner))
            stop("A history-conditional event_learner is required for non-conservative estimation.")
        if (!is.numeric(integration_grid_size) || length(integration_grid_size) != 1L ||
            !is.finite(integration_grid_size) || integration_grid_size < 2 ||
            integration_grid_size != floor(integration_grid_size))
            stop("integration_grid_size must be an integer of at least 2.")
    }
    if (length(list(...))) stop("Unused continuous-time arguments: ", paste(names(list(...)), collapse = ", "))
    if (length(subsets)) stop("Continuous-time subset analyses are not yet supported; subset the input data first.")
    if (is.null(x$continuous_prepared) || is.null(x$models$continuous))
        stop("Run prepare_rtmle_data and model_formula before run_rtmle.")
    if (!identical(x$continuous_baseline_snapshot, x$data$baseline_data))
        stop("Baseline data changed since preparation; rerun prepare_rtmle_data and model_formula.")
    if (is.null(time_horizon)) time_horizon <- max(x$time_grid)
    if (!is.numeric(time_horizon) || !length(time_horizon) || any(!is.finite(time_horizon)) ||
        anyDuplicated(time_horizon) || !all(time_horizon %in% x$time_grid[-1L]))
        stop("time_horizon must select positive time-grid node indices.")
    # Public time horizons are grid-node indices; event histories use time values.
    horizon_times <- x$time_grid_scale[time_horizon + 1L]
    if (!length(x$targets)) stop("Define at least one target with target().")
    if (is.null(targets)) targets <- names(x$targets)
    if (!is.character(targets) || !length(targets) || !all(targets %in% names(x$targets))) stop("Unknown targets.")
    regimes <- unique(unlist(lapply(x$targets[targets], `[[`, "regimes")))
    if (!all(regimes %in% names(x$regimes))) stop("A requested regime is not defined.")
    if (any(vapply(x$regimes[regimes], function(regime) !isTRUE(regime$dynamic_intervention), logical(1))))
        continuous_check_static_support(x$continuous_treatment_support, horizon_times)
    if (!is.null(intervention_hook) && !is.function(intervention_hook)) stop("intervention_hook must be a function.")
    for (flag in list(refit, keep_influence, save_fitted_objects,
                      stratify_sequential_regression, intervene_all_sequential_regression))
        if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) stop("Analysis flags must be TRUE or FALSE.")
    if (intervene_all_sequential_regression &&
        any(vapply(x$regimes[regimes], function(regime) isTRUE(regime$dynamic_intervention), logical(1))))
        stop("intervene_all_sequential_regression is not supported for dynamic continuous-time regimes.")
    outcome_learner <- if (is.null(learner)) "learn_oipcw" else learner
    parsed <- continuous_nuisance_learner(outcome_learner)
    parsed_grid <- if (!conservative && !missing(grid_learner)) continuous_nuisance_learner(grid_learner) else parsed
    for (spec in list(parsed, parsed_grid)) {
        if (any(vapply(spec$learners, is_ipcw_glm_learner, logical(1)))) {
            stop("ipcw_glm must be specified as a single continuous-time outcome or grid learner.")
        }
    }
    weighted_outcome <- is_ipcw_glm_learner(parsed)
    fit_propensity <- estimator != "sequential_regression"
    if (fit_propensity && weighted_outcome && missing(propensity_learner)) propensity_learner <- "learn_glm"
    parsed_propensity <- if (fit_propensity) continuous_nuisance_learner(propensity_learner)
    if (is_ipcw_glm_learner(parsed_propensity)) stop("ipcw_glm is an outcome learner.")
    prepared <- data.table::copy(x$continuous_prepared)
    model_options <- x$models$continuous
    # Cache only nuisance fits, which do not depend on the outcome learner.
    specification <- list(propensity = parsed_propensity,
                           censoring = if (!is.null(censoring_learner)) continuous_censoring_specification(censoring_learner),
                           lag = lag, gbound = if (fit_propensity) gbound, model_options = model_options,
                           exclude_latest_covariate = if (fit_propensity) exclude_latest_covariate,
                           seed = seed, save_fitted_objects = save_fitted_objects)
    if (!conservative) specification$event_learner <- continuous_nuisance_learner(event_learner)
    if (weighted_outcome || !conservative) specification$keep_censoring_prediction_models <- TRUE
    # Changing the outcome analysis invalidates results, but can leave the
    # nuisance fits reusable. Changing their specification invalidates both.
    analysis_specification <- list(nuisance = specification, learner = parsed,
                                    estimator = estimator,
                                    stratify = stratify_sequential_regression,
                                    intervene_all = intervene_all_sequential_regression,
                                    intervention_hook = intervention_hook)
    if (!conservative) {
        analysis_specification$integration_grid_size <- integration_grid_size
        analysis_specification$grid_learner <- parsed_grid
    }
    has_censoring <- any(prepared$info$is_censored)
    experimental <- continuous_experimental_features(analysis_specification, x$regimes[regimes], has_censoring)
    if (length(experimental)) {
        message("The continuous-time RTMLE implementation is highly experimental for ",
                paste(experimental, collapse = ", "),
                ". These features lack a successful comparison with contICEIPCW.")
    }
    if (!identical(analysis_specification, x$continuous_analysis_specification)) {
        x$estimate <- x$IC <- x$sequential_outcome_regression <- NULL
        x$run_time_horizons <- NULL
    }
    if (refit || !identical(specification, x$continuous_nuisance_specification)) x$continuous_nuisance <- NULL
    results <- list()
    diagnostics <- list()
    for (regime_name in regimes) {
        regime <- x$regimes[[regime_name]]
        regime$continuous_baseline_covariates <- x$continuous_prepared$baseline_covariates
        regime$continuous_time_covariates <- x$continuous_prepared$time_covariates
        a <- regime$static_intervention
        if (is.null(x$continuous_nuisance[[regime_name]])) {
            x$continuous_nuisance[[regime_name]] <- continuous_propensity(
                prepared, learner = if (fit_propensity) propensity_learner,
                censoring_learner = censoring_learner, lag = lag, static_intervention = a,
                exclude_latest_covariate = specification$exclude_latest_covariate,
                gbound = specification$gbound, seed = seed, save_fitted_objects = save_fitted_objects,
                regime = regime, model_options = model_options,
                event_learner = if (!conservative) event_learner,
                keep_censoring_prediction_models = weighted_outcome)
            if (!save_fitted_objects) {
                x$continuous_nuisance[[regime_name]]$marginal_censoring_fit <- NULL
                x$continuous_nuisance[[regime_name]]$models$censoring <- lapply(
                    x$continuous_nuisance[[regime_name]]$models$censoring,
                    function(model) { model$fit <- NULL; model })
            }
        }
        nuisance <- x$continuous_nuisance[[regime_name]]
        diagnostics <- continuous_merge_diagnostics(diagnostics, nuisance$diagnostics)
        for (model in nuisance$models$treatment)
            diagnostics <- continuous_merge_diagnostics(diagnostics, model$diagnostics)
        # Nuisance models are shared across horizons for the same regime.
        for (node in time_horizon) {
            th <- x$time_grid_scale[node + 1L]
            grid <- if (!conservative && has_censoring)
                seq(0, th, length.out = integration_grid_size + 1)
            if (verbose || progressbar) message("Continuous time: ", regime_name, ", horizon ", th)
            fitted <- continuous_sequential_regression(
                nuisance, th, a, parsed, estimator, lag, seed,
                x$tuning_parameters$prediction_range,
                stratify_sequential_regression, intervene_all_sequential_regression,
                intervention_hook, save_fitted_objects, regime,
                model_options = model_options, q_grid = grid,
                grid_learner = parsed_grid)
            if (!is.null(grid)) {
                phi_c <- withr::with_preserve_seed(
                    continuous_censoring_correction(nuisance, th, grid, fitted$q_grid, fitted$weights))
                fitted$ic <- fitted$ic + phi_c
                fitted$estimate <- fitted$ice_ipcw + mean(fitted$ic)
                fitted$se <- stats::sd(fitted$ic) / sqrt(length(fitted$ic))
            }
            diagnostics <- continuous_merge_diagnostics(diagnostics, fitted$diagnostics)
            for (target_name in targets) {
                if (!regime_name %in% x$targets[[target_name]]$regimes) next
                results[[length(results) + 1L]] <- data.table::data.table(
                    Target = factor(target_name), Regime = factor(regime_name),
                    Target_parameter = factor("Risk"), Time_horizon = node,
                    Estimator = factor(estimator), Estimate = fitted$estimate,
                    P_value = NA_real_, Standard_error = fitted$se,
                    Lower = pmax(0, fitted$estimate - stats::qnorm(.975) * fitted$se),
                    Upper = pmin(1, fitted$estimate + stats::qnorm(.975) * fitted$se),
                    ICE_IPCW = fitted$ice_ipcw, IPW = fitted$ipw)
                horizon_name <- paste0("time_horizon_", node)
                # summary.rtmle uses these subject-aligned curves for contrasts.
                x$IC[[target_name]][[regime_name]][[horizon_name]] <- fitted$ic
                x$sequential_outcome_regression[[target_name]][[regime_name]][[horizon_name]] <- fitted$models
            }
        }
    }
    # Replace only the requested result rows when adding horizons to a fit.
    new_results <- data.table::rbindlist(results)
    old <- x$estimate$Main_analysis
    if (!is.null(old)) {
        keys <- c("Target", "Regime", "Time_horizon", "Estimator")
        old <- old[!new_results, on = keys]
        new_results <- data.table::rbindlist(list(old, new_results), fill = TRUE)
    }
    data.table::setkeyv(new_results, c("Target", "Regime", "Target_parameter", "Time_horizon", "Estimator"))
    x$estimate$Main_analysis <- new_results
    x$learner <- parsed
    x$unparsed_learner <- outcome_learner
    x$continuous_nuisance_specification <- specification
    x$continuous_analysis_specification <- analysis_specification
    x$run_time_horizons <- sort(unique(c(x$run_time_horizons, time_horizon)))
    x$save_fitted_objects <- save_fitted_objects
    x$diagnostics <- diagnostics
    x
}
