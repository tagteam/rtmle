library(testthat)
library(rtmle)

test_that("changing the learner clears previous analysis state", {
    data(rtmle_object)
    x <- data.table::copy(rtmle_object)
    x$estimate$Old_analysis <- data.table::copy(x$estimate$Main_analysis)

    x <- run_rtmle(x,
                   learner = "learn_glmnet",
                   time_horizon = 1,
                   seed = 1,
                   verbose = FALSE)

    expect_equal(x$learner$name, "learn_glmnet")
    expect_equal(sort(unique(x$estimate$Main_analysis$Time_horizon)), 1)
    expect_equal(names(x$estimate), "Main_analysis")
    expect_equal(x$run_time_horizons, 1)
    expect_length(x$regimes$Always_A$ipw_last_nodes, 1)

    previous_model <- x$models$time_1$outcome[[paste0(x$names$outcome, "_2")]]
    expect_false(any(c("fit", "fit_summary") %in% names(previous_model)))
})
