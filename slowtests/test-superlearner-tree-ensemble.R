library(testthat)
library(rtmle)

tree_learners <- list(
    forest = list(fun = "learn_ranger", num.trees = 3L, replace = FALSE,
                  sample.fraction = 1, min.node.size = 1L, num.threads = 1L),
    boost = list(fun = "learn_xgboost", xgboost.max_depth = 1L,
                 xgboost.nrounds = 1L, xgboost.learning_rate = 1,
                 xgboost.reg_lambda = 0, xgboost.nthreads = 1L)
)

test_that("a cross-validated ranger and xgboost ensemble recovers known group means", {
    skip_if_not_installed("ranger")
    skip_if_not_installed("xgboost", minimum_version = "3.1.2")

    d <- data.table::data.table(id = 1:80, x = rep(0:1, each = 40),
                                y = c(rep(c(0, .25, .25, .5), 10),
                                      rep(c(.5, .75, .75, 1), 10)))
    # A full-sample forest and an unregularized single boosting step each
    # predict the training mean conditional on x in every cross-validation fold.
    parsed <- parse_learners(list(folds = 3, ensemble_method = "ipa", learners = tree_learners))
    fit <- superlearn(folds = 3, seed = 13, learners = parsed$learners,
                      parse_learners = FALSE, character_formula = "y ~ x", outcome_variable = "y",
                      id_variable = "id", data = d, intervened_data = data.table::data.table(x = c(1, 0, 0)),
                      ensemble_method = "ipa", save_fitted_objects = TRUE)
    expect_equal(fit$fit_summary, c(forest = .5, boost = .5), tolerance = 1e-6)
    expect_equal(as.numeric(fit$predicted_values), c(.75, .25, .25), tolerance = 1e-6)
    expect_s3_class(fit$fit$forest, "ranger")
    expect_s3_class(fit$fit$boost, "xgb.Booster")
})

test_that("TMLE with a tree ensemble recovers both one-interval risks and IC standard errors", {
    skip_if_not_installed("ranger")
    skip_if_not_installed("xgboost", minimum_version = "3.1.2")
    id <- seq_len(80)
    x <- rtmle_init(time_grid = 0:1, name_id = "id", name_outcome = "Y")
    x <- add_baseline_data(x, data.frame(id = id))
    treatment <- data.frame(id = id, A_0 = rep(0:1, each = 40))
    outcome <- data.frame(id = id, Y_1 = c(rep(c(0, 0, 0, 1), 10), rep(c(0, 1, 1, 1), 10)))
    x <- add_wide_data(x, outcome_data = outcome, timevar_data = list(A = treatment))
    x <- prepare_rtmle_data(x, verbose = FALSE)
    x <- regime(x, name = "treated", treatment_variables = "A", intervention = 1, verbose = FALSE)
    x <- regime(x, name = "untreated", treatment_variables = "A", intervention = 0, verbose = FALSE)
    x <- target(x, name = "risk", regimes = c("treated", "untreated"))
    x <- model_formula(x, verbose = FALSE)
    fit <- run_rtmle(x, time_horizon = 1L, seed = 13, save_fitted_objects = TRUE,
                     learner = list(folds = 3, ensemble_method = "ipa", learners = tree_learners))

    for (regime in c("treated", "untreated")) {
        outcome_fit <- fit$models$time_0$outcome$Y_1$fit[[regime]]$sequence_time_1
        expect_s3_class(outcome_fit$fit$forest, "ranger")
        expect_s3_class(outcome_fit$fit$boost, "xgb.Booster")
        expect_equal(outcome_fit$fit_summary, c(forest = .5, boost = .5), tolerance = 1e-6)
    }

    risks <- fit$estimate$Main_analysis
    rows <- match(c("treated", "untreated"), as.character(risks$Regime))
    expect_equal(risks$Estimate[rows], c(30, 10) / 40, tolerance = 1e-6)
    # In both arms, IC = 2 * I(A = a) * (Y - risk), with squared sum
    # 4 * 40 * (3/4) * (1/4) = 30. The empirical variance uses n - 1.
    expect_equal(risks$Standard_error[rows], rep(sqrt(30 / (80 * 79)), 2), tolerance = 1e-6)
})
