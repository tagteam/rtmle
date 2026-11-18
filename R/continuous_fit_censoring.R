# Fit one absolute-time hazard model, shared by all event indices and horizons.
continuous_fit_censoring <- function(prepared, specification, model_options, call_learner) {
    training <- continuous_censoring_intervals(prepared)
    predictors <- unique(c(prepared$baseline_covariates, prepared$time_covariates))
    response <- "survival::Surv(.rtmle_start, .rtmle_stop, .rtmle_censored)"
    custom <- specification$current_formula
    if (is.null(custom)) {
        # Use current values, alongside fixed baseline predictors. Their initial
        # *_0 versions remain available when explicitly requested in a formula.
        terms <- c(setdiff(prepared$baseline_covariates, paste0(prepared$time_covariates, "_0")),
                   prepared$time_covariates)
        formula <- continuous_nuisance_formula(response, terms)
        intercept <- 1L
    } else {
        formula <- paste(response, custom)
        # Resolve '.' using only model columns, excluding the subject identifier.
        model_data <- training[, c(".rtmle_start", ".rtmle_stop", ".rtmle_censored", predictors), with = FALSE]
        parsed <- stats::terms(stats::as.formula(formula), data = model_data, keep.order = TRUE)
        terms <- attr(parsed, "term.labels")
        offsets <- attr(parsed, "offset")
        if (length(offsets)) terms <- c(terms, vapply(offsets, function(i) {
            paste(deparse(attr(parsed, "variables")[[i + 1L]]), collapse = " ")
        }, character(1)))
        intercept <- attr(parsed, "intercept")
        used <- all.vars(stats::delete.response(parsed))
        if (!all(used %in% predictors))
            stop("censoring_formula must use baseline or current time-varying covariates.")
    }
    name <- model_options$node_names$censoring
    if (is.null(name)) name <- "C"
    choices <- lapply(seq_len(max(prepared$info$last_event)), function(k) {
        continuous_apply_model_rules(terms, union(predictors, terms), names(training),
                                     model_options, paste0(name, "_", k), prepared$time_covariates)
    })
    if (any(!vapply(choices, function(x) setequal(x, choices[[1L]]), logical(1))))
        stop("A time-varying censoring model is shared across event indices; censoring model rules must select the same predictors.")
    terms <- choices[[1L]]
    if (is.null(custom)) {
        formula <- continuous_nuisance_formula(response, terms)
    } else {
        if (!intercept) terms <- c("0", terms)
        # Terms can include prop(), const(), strata(), interactions and offsets.
        rhs <- if (length(terms)) paste(terms, collapse = " + ") else "1"
        formula <- paste(response, "~", rhs)
    }
    fitted <- call_learner(c(list(character_formula = formula, data = training,
        intervened_data = training, time_variable = ".rtmle_stop", time_ref = ".rtmle_start",
        save_fitted_objects = TRUE), specification$learner$args))
    continuous_nuisance_predictions(fitted, nrow(training), "Censoring", TRUE)
    list(fitted = fitted, intervals = training, formula = formula)
}
