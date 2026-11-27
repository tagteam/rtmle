# Keep this list aligned with the contICEIPCW comparisons in the sandbox and
# tests/testthat/test-continuous-source-scenarios.R.
continuous_experimental_features <- function(analysis, regimes, has_censoring) {
    is_learner <- function(specification, name, fun) {
        identical(specification$fun, name) || identical(specification$fun, fun)
    }

    features <- character()
    if (analysis$estimator == "tmle") features <- c(features, "continuous-time TMLE")
    if (has_censoring && !is.null(analysis$integration_grid_size))
        features <- c(features, "the non-conservative censoring correction")
    if (any(vapply(regimes, function(regime) isTRUE(regime$dynamic_intervention), logical(1))))
        features <- c(features, "dynamic treatment regimes")
    if (!is.null(analysis$intervention_hook)) features <- c(features, "custom intervention hooks")

    outcome <- analysis$learner
    link <- outcome$args[["link"]]
    lambda <- outcome$args[["lambda"]]
    gaussian <- is_learner(outcome, "learn_glmnet", learn_glmnet) &&
        identical(outcome$args[["family"]], "gaussian") && is.numeric(lambda) && isTRUE(lambda == 0)
    if (!is_learner(outcome, "learn_oipcw", learn_oipcw) && !gaussian &&
        !(is_ipcw_glm_learner(outcome) && (is.null(link) || isTRUE(link %in% c("logit", "probit")))))
        features <- c(features, paste0("outcome learner '", outcome$name, "'"))

    nuisance <- analysis$nuisance
    propensity <- nuisance$propensity
    if (!is.null(propensity)) {
        if (!is.null(nuisance$gbound)) features <- c(features, "propensity bounding (gbound)")
        if (any(lengths(lapply(regimes, `[[`, "adherence_model_strata")) > 0L))
            features <- c(features, "adherence model strata")

        family <- propensity$args[["family"]]
        logistic_glmnet <- is_learner(propensity, "learn_glmnet", learn_glmnet) &&
            identical(propensity$args[["selector"]], "min") && is.null(propensity$args[["lambda"]]) &&
            (is.null(family) || identical(family, "binomial"))
        if (!is_learner(propensity, "learn_glm", learn_glm) && !logistic_glmnet)
            features <- c(features, paste0("propensity learner '", propensity$name, "'"))
    }

    censoring <- nuisance$censoring
    if (has_censoring && !is.null(censoring)) {
        if (censoring$time_varying) features <- c(features, "shared time-varying censoring models")
        if (!is_learner(censoring$learner, "learn_coxph", learn_coxph)) {
            features <- c(features, paste0("censoring learner '", censoring$learner$name, "'"))
        } else if (!is.null(censoring$formula) && grepl("\\bstrata[[:space:]]*\\(", censoring$formula)) {
            features <- c(features, "stratified marginal Cox censoring")
        }
    }

    formula_options <- c("Markov", "exclude_variables", "inclusion_rules", "exclusion_rules")
    if (any(lengths(nuisance$model_options[formula_options]) > 0L))
        features <- c(features, "custom model_formula() rules")
    features
}
