test_that("fits covered by the contICEIPCW comparisons are quiet by default", {
    x <- continuous_analysis()
    expect_silent(run_rtmle(x, estimator = "one-step"))
    expect_silent(run_rtmle(x, estimator = "one-step", conservative = FALSE,
                           integration_grid_size = 2L))
    expect_silent(run_rtmle(x, estimator = "sequential_regression", gbound = 0.6,
                           learner = list(name = "probit", fun = learn_ipcw_glm, link = "probit")))
})

test_that("one message lists experimental features across horizons and regimes", {
    x <- continuous_analysis(horizons = c(1, 2))
    messages <- capture_messages(run_rtmle(x, estimator = "tmle", gbound = 0.6,
                                           time_horizon = 1:2, verbose = FALSE))
    expect_length(messages, 1L)
    expect_match(messages[[1L]], "highly experimental.*continuous-time TMLE.*propensity bounding.*contICEIPCW")
})

test_that("the experimental message covers censoring features in use", {
    data <- continuous_example()
    data$timevarying_data$event[1L] <- "C"
    data$timevarying_data$time[1L] <- 1
    x <- continuous_analysis(data, regimes = "Always_A")
    censoring <- list(name = "shared", fun = "learn_coxph", time_varying = TRUE,
                      censoring_formula = ~ 1)
    expect_message(run_rtmle(x, estimator = "one-step", conservative = FALSE,
                            integration_grid_size = 2L, censoring_learner = censoring),
                   "highly experimental.*non-conservative censoring correction.*shared time-varying censoring")
})
