continuous_learner_analysis <- function() {
    set.seed(34)
    data <- simulate_continuous_time_data(300, uncensored = FALSE,
                                           static_intervention_baseline = NULL)
    continuous_analysis(data, 720,
                        regimes = "Always_A", min_events = 10)
}

test_that("continuous outcome and propensity models accept glmnet specifications", {
    learner <- list(name = "penalized", fun = "learn_glmnet", lambda = 0.05)
    fit <- run_rtmle(continuous_learner_analysis(), learner = learner,
                     estimator = "one-step", seed = 31,
                     save_fitted_objects = TRUE)
    result <- fit$estimate$Main_analysis
    expect_true(all(is.finite(result$Estimate)))
    expect_true(all(is.finite(result$Standard_error)))
    expect_s3_class(fit$continuous_nuisance$Always_A$models$treatment$A_0$fit,
                    "glmnet")
    outcome_models <- fit$sequential_outcome_regression$risk$Always_A$time_horizon_1
    expect_true(length(outcome_models) > 1L)
    for (model in outcome_models) expect_s3_class(model$fit, "glmnet")
})

test_that("continuous models accept forest specifications and propensity bounds", {
    skip_if_not_installed("ranger")
    learner <- list(name = "forest", fun = "learn_ranger", num.trees = 20,
                    min.node.size = 5, num.threads = 1, seed = 31)
    fit <- run_rtmle(continuous_learner_analysis(), learner = learner,
                     estimator = "one-step", gbound = 0.05,
                     save_fitted_objects = TRUE)
    expect_true(all(is.finite(fit$estimate$Main_analysis$Estimate)))
    propensity <- fit$continuous_nuisance$Always_A
    expect_s3_class(propensity$models$treatment$A_0$fit, "ranger")
    columns <- grep("^propensity_", names(propensity$data), value = TRUE)
    for (column in columns) {
        probabilities <- stats::na.omit(propensity$data[[column]])
        expect_true(all(probabilities >= 0.05 & probabilities <= 1))
    }
    outcome_models <- fit$sequential_outcome_regression$risk$Always_A$time_horizon_1
    for (model in outcome_models) expect_s3_class(model$fit, "ranger")
})

test_that("continuous models reuse the discrete pipeline's ensemble interface", {
    learner <- list(folds = 2, ensemble_method = "discrete", learners = list(
        list(name = "penalized", fun = "learn_glmnet", lambda = 0.05),
        list(name = "more_penalized", fun = "learn_glmnet", lambda = 0.1)
    ))
    x <- continuous_learner_analysis()
    fit <- run_rtmle(x, learner = learner, estimator = "one-step", seed = 31,
                     save_fitted_objects = TRUE)
    expect_identical(fit$learner$name, "superlearn")
    expect_equal(length(fit$learner$learners), 2L)
    expect_true(all(is.finite(fit$estimate$Main_analysis$Estimate)))
    expect_true(all(is.finite(fit$continuous_nuisance$Always_A$data$propensity_0)))
    outcome_models <- fit$sequential_outcome_regression$risk$Always_A$time_horizon_1
    for (model in outcome_models) {
        expect_true(length(model$fit) > 0L)
        for (member in model$fit) expect_s3_class(member, "glmnet")
    }
    baseline_members <- fit$continuous_nuisance$Always_A$models$treatment$A_0$fit
    expect_true(length(baseline_members) > 0L)
    for (member in baseline_members) expect_s3_class(member, "glmnet")
    repeated <- run_rtmle(x, learner = learner, estimator = "one-step", seed = 31,
                          save_fitted_objects = TRUE)
    expect_equal(repeated$estimate$Main_analysis, fit$estimate$Main_analysis)
})
