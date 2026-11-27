test_that("the weighted learner uses ordinary quasi-binomial GLM links", {
    d <- data.frame(Y = c(0, .2, .8, 1), A = c(0, 0, 1, 1), w = c(0, 2, 1, 3))
    for (link in c("logit", "probit")) {
        fit <- learn_ipcw_glm("Y ~ A", d, d, weights = d$w, link = link,
                              save_fitted_objects = TRUE)
        expected <- stats::glm(Y ~ A, data = d, weights = w, family = quasibinomial(link))
        expect_equal(fit$predicted_values, unname(stats::fitted(expected)))
        expect_identical(fit$fit$family$link, link)
        reused <- learn_ipcw_glm("Y ~ A", data = NULL, intervened_data = d[4, ], reuse_fit = fit)
        expect_equal(reused$predicted_values, fit$predicted_values[4])
    }
})

test_that("IPCW GLMs use censoring weights at the horizon and at the grid cutoff", {
    d <- list(baseline_data = data.frame(id = 1:6, A_0 = 1, L_0 = 0),
              timevarying_data = data.frame(id = 1:6, time = c(.75, 1, 1, 2, 3, 3),
                  event = c("D", "C", "Y", "Y", "tauend", "tauend"), A = 1, L = 0))
    x <- continuous_analysis(d, c(1, 2), regimes = "Always_A", min_events = 0)

    # The Cox censoring hazard jumps by 1/5 at time 1. At that cutoff the
    # tied outcome uses the left limit. The three later events contribute
    # zeros weighted by the inverse right limit, and dropout has zero weight.
    b <- exp(1/5)
    q <- c(1, 1 + b)/(2 + 3*b)
    for (censoring in list(
        list(name = "marginal", fun = "learn_coxph", marginal_censoring_formula = ~ 1),
        list(name = "conditional", fun = "learn_coxph", marginal_censoring = FALSE),
        list(name = "shared", fun = "learn_coxph", time_varying = TRUE, censoring_formula = ~ 1)
    )) {
        conservative <- run_rtmle(x, learner = "ipcw_glm", estimator = "one-step",
                                  censoring_learner = censoring, time_horizon = 1:2,
                                  save_fitted_objects = TRUE)
        expect_equal(conservative$estimate$Main_analysis$ICE_IPCW, q)
        expect_identical(conservative$continuous_nuisance_specification$propensity$fun, "learn_glm")
        glm <- conservative$sequential_outcome_regression$risk$Always_A$time_horizon_1$event_1$fit
        expect_equal(unname(glm$prior.weights), c(1, 0, 1, b, b, b))

        updated <- run_rtmle(conservative, learner = "ipcw_glm", estimator = "one-step",
                            censoring_learner = censoring, time_horizon = 2, conservative = FALSE,
                            integration_grid_size = 2L, refit = FALSE)
        difference <- updated$IC$risk$Always_A$time_horizon_2 -
            conservative$IC$risk$Always_A$time_horizon_2
        # Event survival just before dropout includes the competing event at .75.
        correction <- diff(q) * (c(0, exp(1/6), 0, 0, 0, 0)*b - c(0, 1, 1, 1, 1, 1)*(b - 1))
        expect_equal(difference, correction)
    }
})
