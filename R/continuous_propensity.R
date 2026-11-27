continuous_propensity <- function(prepared_data, learner = "learn_glm",
                                  censoring_learner = "learn_coxph", lag = NULL,
                                  static_intervention = 1,
                                  exclude_latest_covariate = NULL,
                                  gbound = NULL, seed = NULL,
                                  save_fitted_objects = FALSE,
                                  regime = NULL, model_options = NULL, event_learner = NULL,
                                  keep_censoring_prediction_models = FALSE) {
    if (!is.list(prepared_data) || is.null(prepared_data$wide_data) ||
        is.null(prepared_data$info)) stop("prepared_data must come from prepare_continuous_data.")
    if (length(static_intervention) != 1L || is.na(static_intervention) ||
        !static_intervention %in% c(0, 1)) stop("static_intervention must be 0 or 1.")
    if (!is.null(gbound) && (!is.numeric(gbound) || length(gbound) != 1L ||
        !is.finite(gbound) || gbound <= 0 || gbound > 1)) stop("gbound must lie in (0, 1].")
    if (!is.null(lag) && (!is.numeric(lag) || length(lag) != 1L ||
        !is.finite(lag) || lag < 1 || lag != floor(lag))) stop("lag must be a positive integer.")
    if (!is.null(regime)) {
        if (is.null(regime$continuous_baseline_covariates))
            regime$continuous_baseline_covariates <- prepared_data$baseline_covariates
        if (is.null(regime$continuous_time_covariates))
            regime$continuous_time_covariates <- prepared_data$time_covariates
    }
    data <- data.table::copy(prepared_data$wide_data)
    info <- prepared_data$info
    treatment_learner <- if (!is.null(learner)) continuous_nuisance_learner(learner)
    has_censoring <- any(info$is_censored)
    marginal_fit <- NULL
    models <- list(treatment = list(), censoring = list())
    learner_warnings <- NULL
    censoring_name <- model_options$node_names$censoring
    if (is.null(censoring_name) || !length(censoring_name)) censoring_name <- "C"
    censoring <- continuous_censoring_specification(censoring_learner)
    event_model <- if (!is.null(event_learner)) continuous_event_specification(event_learner)
    needs_grid <- !is.null(event_model) && has_censoring
    needs_cutoff <- has_censoring && (needs_grid || keep_censoring_prediction_models)
    grid <- if (needs_grid) list(event_learner = event_model)
    cutoff_models <- if (needs_cutoff) list(specification = censoring)
    if (has_censoring) {
        if (is.null(censoring)) stop("A censoring learner is required when censoring is present.")
        hazard_learner <- censoring$learner
    }

    call_hazard <- function(args, specification = hazard_learner) {
        withCallingHandlers(do.call(specification$fun, args), warning = function(w) {
            learner_warnings <<- rbind(learner_warnings,
                                      data.table::data.table(learner = specification$name,
                                                             warning = conditionMessage(w)))
            invokeRestart("muffleWarning")
        })
    }
    shared_censoring <- NULL
    if (has_censoring && censoring$time_varying) {
        shared_censoring <- continuous_fit_censoring(prepared_data, censoring, model_options, call_hazard)
        models$censoring$time_varying <- shared_censoring$fitted
        if (needs_cutoff) cutoff_models$shared <- shared_censoring
    }
    # A marginal model is fitted once using baseline predictors. Its saved fit
    # is then used to predict survival across each person's event intervals.
    # Conditional models, in contrast, are fitted separately at each interval.
    if (has_censoring && censoring$marginal) {
        training <- prepared_data$data_marginal_censoring
        censoring_nodes <- paste0(censoring_name, "_", seq_len(max(info$last_event)))
        if (!is.null(censoring$formula)) {
            # The optional formula can contain terms such as strata(group),
            # which must remain formula terms rather than quoted column names.
            custom_formula <- paste('survival::Surv(time, event == "C")', censoring$formula)
            # Resolve `.` against the survival formula and its training data;
            # terms() stores offsets separately from its ordinary term labels.
            rhs_terms <- stats::terms(stats::as.formula(custom_formula),
                                      data = training, keep.order = TRUE)
            offsets <- attr(rhs_terms, "offset")
            variables <- attr(rhs_terms, "variables")
            offset_terms <- vapply(offsets, function(i) {
                paste(deparse(variables[[i + 1L]]), collapse = " ")
            }, character(1))
            rhs <- c(attr(rhs_terms, "term.labels"), offset_terms)
            intercept <- attr(rhs_terms, "intercept")
            available <- union(prepared_data$baseline_covariates, rhs)
            all_names <- union(names(training), rhs)
            by_node <- lapply(censoring_nodes, function(node) {
                continuous_apply_model_rules(rhs, available, all_names,
                                             model_options, node, prepared_data$time_covariates)
            })
        } else {
            by_node <- lapply(censoring_nodes, function(node) {
                continuous_history(training, prepared_data$time_covariates,
                                   prepared_data$baseline_covariates,
                                   "marginal", k = 1L, model_options = model_options,
                                   outcome_name = node)
            })
        }
        # A marginal model is fitted once, so it cannot obey conflicting
        # specifications for different censoring nodes.
        if (any(!vapply(by_node, function(terms) setequal(terms, by_node[[1L]]), logical(1))))
            stop("Censoring rules select different predictors across event indices, but a marginal censoring model is shared. Use consistent rules or a conditional censoring model.")
        covariates <- by_node[[1L]]
        if (!is.null(censoring$formula)) {
            if (setequal(covariates, rhs)) {
                # Preserve the user's exact formula, including its dot and
                # offset expressions, if no model rule changed the terms.
                formula <- custom_formula
            } else {
                if (!intercept) covariates <- c("0", covariates)
                rhs <- if (length(covariates)) paste(covariates, collapse = " + ") else "1"
                formula <- paste('survival::Surv(time, event == "C") ~', rhs)
            }
        } else formula <- continuous_nuisance_formula('survival::Surv(time, event == "C")', covariates)
        marginal_fit <- call_hazard(c(list(character_formula = formula, data = training,
                                           intervened_data = training, time_variable = "time",
                                           time_ref = NULL, save_fitted_objects = TRUE),
                                      hazard_learner$args))
        continuous_nuisance_predictions(marginal_fit, nrow(training), "Censoring", TRUE)
        marginal_formula <- formula
        models$censoring$marginal <- marginal_fit
        if (needs_cutoff) cutoff_models$marginal <- list(fit = marginal_fit, formula = formula)
    }

    # At each event index, predict the probability of remaining uncensored
    # from the previous event to the next. The extra pooled endpoint is used
    # when the analysis combines later events at a horizon.
    last_event <- max(info$last_event)
    for (k in seq_len(last_event)) {
        previous_event <- data[[paste0("event_", k - 1L)]]
        rows <- if (k == 1L) seq_len(nrow(data)) else which(previous_event %in% c("A", "L"))
        survival_columns <- paste0("survival_censoring_", k)
        pooled <- k %in% info$last_event
        if (pooled) survival_columns <- c(survival_columns, paste0("survival_censoring_pooled_", k))
        for (column in survival_columns) data.table::set(data, j = column, value = rep(NA_real_, nrow(data)))
        if (!length(rows)) next
        for (column in survival_columns) {
            pooled_column <- startsWith(column, "survival_censoring_pooled_")
            fit_grid_model <- needs_grid && (k < last_event || pooled_column)
            endpoint <- paste0(if (pooled_column) "time_pooled_" else "time_", k)
            event <- paste0(if (pooled_column) "event_pooled_" else "event_", k)
            prediction <- data.table::copy(data[rows])
            if (!has_censoring) {
                probabilities <- rep(1, length(rows))
            } else if (censoring$time_varying) {
                probabilities <- continuous_predict_censoring(shared_censoring, prediction, endpoint,
                    paste0("time_", k - 1L), censoring, call_hazard)
            } else if (censoring$marginal) {
                # The learner must reuse the marginal fit and predict survival
                # over the requested interval for each prediction row.
                predicted <- call_hazard(c(list(character_formula = marginal_formula,
                    data = training, intervened_data = prediction, time_variable = endpoint,
                    time_ref = paste0("time_", k - 1L), save_fitted_objects = TRUE,
                    reuse_fit = marginal_fit), hazard_learner$args))
                probabilities <- continuous_nuisance_predictions(predicted, length(rows), "Censoring", TRUE)
            } else {
                prediction <- continuous_gap_history(prediction, k, endpoint)
                covariates <- continuous_history(prediction, prepared_data$time_covariates,
                                                  prepared_data$baseline_covariates, "hazard", lag, k,
                                                  model_options = model_options,
                                                  outcome_name = paste0(censoring_name, "_", k))
                outcome <- paste0('survival::Surv(`', endpoint, '`, `', event, '` == "C")')
                formula <- continuous_nuisance_formula(outcome, covariates)
                fitted <- call_hazard(c(list(character_formula = formula, data = prediction,
                    intervened_data = prediction, time_variable = endpoint, time_ref = NULL,
                    save_fitted_objects = TRUE), hazard_learner$args))
                probabilities <- continuous_nuisance_predictions(fitted, nrow(prediction), "Censoring", TRUE)
                models$censoring[[column]] <- fitted
                if (needs_cutoff && (k < last_event || pooled_column))
                    cutoff_models$conditional[[endpoint]] <- list(formula = formula, fit = fitted,
                                                                   training = prediction, rows = rows)
            }
            if (fit_grid_model) {
                # Model the gap duration given the history at entry, including
                # the original timing of earlier events.
                if (!censoring$marginal && !censoring$time_varying) prediction <- data[rows]
                start <- prediction[[paste0("time_", k - 1L)]]
                data.table::set(prediction, j = endpoint, value = prediction[[endpoint]] - start)
                covariates <- continuous_history(prediction, prepared_data$time_covariates,
                    prepared_data$baseline_covariates, "event_survival", lag, k,
                    model_options = model_options, outcome_name = paste0("S_", k))
                outcome <- paste0('survival::Surv(`', endpoint, '`, `', event,
                                  '` %in% c("A", "L", "Y", "D"))')
                formula <- continuous_nuisance_formula(outcome, covariates)
                # Adding S_k must not change the random draws used by the
                # existing treatment, censoring and horizon regressions.
                fitted <- withr::with_preserve_seed(call_hazard(c(list(
                    character_formula = formula, data = prediction,
                    intervened_data = prediction, time_variable = endpoint, time_ref = NULL,
                    save_fitted_objects = TRUE), event_model$args), event_model))
                continuous_nuisance_predictions(fitted, nrow(prediction), "Event survival")
                grid$event[[endpoint]] <- list(formula = formula, fit = fitted,
                                               training = prediction, rows = rows)
            }
            data.table::set(data, i = rows, j = column, value = as.numeric(probabilities))
        }
    }

    # Fit the probability of the observed assignment agreeing with the regime.
    # A natural-course decision always agrees, so its probability is one and
    # those rows do not need an adherence model.
    fit_treatment <- function(rows, outcome, covariates, label) {
        training <- data.table::copy(data[rows])
        strata <- continuous_adherence_strata(regime$adherence_model_strata, label, names(training),
                                               prepared_data$baseline_covariates,
                                               prepared_data$time_covariates)
        policy <- if (is.null(regime)) list(value = rep(static_intervention, nrow(training)),
                                             natural = rep(FALSE, nrow(training))) else
            continuous_intervention_policy(training, regime, label)
        active <- which(!policy$natural)
        probabilities <- rep(1, nrow(training))
        if (!length(active)) return(list(probabilities = probabilities, fitted = NULL))
        data.table::set(training, j = "rtmle_treatment_outcome",
                        value = as.numeric(training[[outcome]] == policy$value))
        if (length(strata) && any(vapply(training[active, strata, with = FALSE], anyNA, logical(1))))
            stop("`adherence_model_strata` variables must be observed for every row used to fit an adherence propensity.")
        groups <- if (length(strata)) {
            groups_data <- training[active, list(rows = .I), by = strata]
            groups_data[, list(indices = list(rows)), by = strata]$indices
        } else list(active)
        fits <- lapply(groups, function(group) {
            group_covariates <- continuous_varying_covariates(training[group], covariates)
            formula <- continuous_nuisance_formula("rtmle_treatment_outcome", group_covariates)
            fit <- fitter(intervention_node = label, learner = treatment_learner,
                          formula = formula, data = training[group], intervened_data = training[group],
                          id_variable = "id", minority_threshold = 0, seed = seed,
                          save_fitted_objects = save_fitted_objects)
            probabilities[group] <<- continuous_nuisance_predictions(fit, length(group), "Treatment propensity")
            fit
        })
        fitted <- if (!length(strata)) fits[[1L]] else list(
            fit = if (save_fitted_objects) lapply(fits, `[[`, "fit"),
            fit_summary = lapply(fits, `[[`, "fit_summary"),
            diagnostics = Reduce(continuous_merge_diagnostics, lapply(fits, `[[`, "diagnostics"), init = list()))
        if (!is.null(gbound)) probabilities <- pmax(probabilities, gbound)
        list(probabilities = probabilities, fitted = fitted)
    }
    # Baseline treatment is modeled for everyone; later treatment models are
    # fitted only at observed treatment opportunities, using pre-decision history.
    if (!is.null(treatment_learner)) {
        baseline_covariates <- setdiff(prepared_data$baseline_covariates, "A_0")
        baseline_covariates <- continuous_history(data, prepared_data$time_covariates,
                                                   baseline_covariates, "baseline", k = 0L,
                                                   model_options = model_options,
                                                   outcome_name = "A_0")
        baseline <- fit_treatment(seq_len(nrow(data)), "A_0", baseline_covariates, 0L)
        data.table::set(data, j = "propensity_0", value = baseline$probabilities)
        models$treatment$A_0 <- baseline$fitted
        for (k in seq_len(last_event - 1L)) {
            rows <- which(data[[paste0("event_", k)]] == "A")
            column <- paste0("propensity_", k)
            data.table::set(data, j = column, value = rep(NA_real_, nrow(data)))
            if (!length(rows)) next
            covariates <- continuous_history(data[rows], prepared_data$time_covariates,
                                              prepared_data$baseline_covariates, "propensity", lag, k,
                                              exclude_latest_covariate,
                                              model_options = model_options,
                                              outcome_name = paste0("A_", k))
            fitted <- fit_treatment(rows, paste0("A_", k), covariates, k)
            data.table::set(data, i = rows, j = column, value = fitted$probabilities)
            models$treatment[[paste0("A_", k)]] <- fitted$fitted
        }
    }
    out <- list(data = data, prepared_data_object = prepared_data,
                marginal_censoring_fit = marginal_fit, models = models,
                diagnostics = list(learner_warnings = learner_warnings),
                static_intervention = static_intervention)
    if (needs_grid) out$grid_survival <- grid
    if (needs_cutoff) out$censoring_prediction <- cutoff_models
    class(out) <- "debiased_prepared"
    out
}
