test_that("replacing continuous data invalidates nuisance fits even without refit", {
    x <- run_rtmle(continuous_analysis(), estimator = "one-step")
    expect_equal(x$estimate$Main_analysis[Regime == "Always_A", Estimate], 0.5,
                 tolerance = 1e-7)

    data <- continuous_example()
    data$timevarying_data$event <- c("Y", "tauend", "tauend", "tauend")
    x <- continuous_add_example_data(x, data)
    x <- prepare_rtmle_data(x)
    x <- model_formula(x, verbose = FALSE)
    x <- run_rtmle(x, estimator = "one-step", refit = FALSE)

    expect_equal(x$estimate$Main_analysis$Estimate, c(0, 0.5), tolerance = 1e-7)
    expect_equal(x$estimate$Main_analysis$IPW, c(0, 0.5))
    expect_equal(as.character(x$continuous_nuisance$Always_A$data$event_1),
                 data$timevarying_data$event)
})

test_that("an explicitly empty outcome stream clears earlier continuous outcomes", {
    x <- run_rtmle(continuous_analysis(), estimator = "one-step")
    expect_equal(x$estimate$Main_analysis[Regime == "Always_A", Estimate], 0.5,
                 tolerance = 1e-7)

    data <- continuous_example()
    data$timevarying_data$event <- "tauend"
    x <- continuous_add_example_data(x, data)
    expect_equal(NROW(x$long_data$outcome_data), 0L)
    x <- model_formula(prepare_rtmle_data(x), verbose = FALSE)
    x <- run_rtmle(x, estimator = "one-step", refit = FALSE)

    expect_equal(x$estimate$Main_analysis$Estimate, c(0, 0))
    expect_equal(x$estimate$Main_analysis$IPW, c(0, 0))
    expect_true(all(x$continuous_nuisance$Always_A$data$event_1 == "tauend"))
})

test_that("baseline data and continuous events can be added in either order", {
    data <- continuous_example()
    initialize <- function() {
        rtmle_init(c(0, 2), name_id = "id", name_time = "time",
                    name_outcome = "Y", name_competing = "D",
                    name_censoring = "C", continuous_time = TRUE)
    }
    add_events <- function(x) {
        continuous_add_example_events(x, data)
    }
    baseline_first <- add_events(add_baseline_data(initialize(), data$baseline_data))
    events_first <- add_baseline_data(add_events(initialize()), data$baseline_data)
    baseline_first <- prepare_rtmle_data(baseline_first)
    events_first <- prepare_rtmle_data(events_first)

    expect_equal(events_first$prepared_data, baseline_first$prepared_data)
    expect_equal(events_first$continuous_prepared, baseline_first$continuous_prepared)
    expect_equal(events_first$prepared_data$A_0, data$baseline_data$A_0)
})

test_that("replacing only baseline data clears and refits continuous analyses", {
    outcome <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = rep(0.2, nrow(intervened_data)))
    }
    propensity <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = stats::plogis(intervened_data$L_0))
    }
    fit <- function(x, refit = TRUE) {
        run_rtmle(x, learner = outcome, propensity_learner = propensity,
                   estimator = "one-step", refit = refit)
    }
    x <- fit(continuous_analysis())
    before <- x$estimate$Main_analysis[Regime == "Always_A", Estimate]
    expect_equal(before, 0.5, tolerance = 1e-7)

    data <- continuous_example()
    data$baseline_data$L_0 <- c(-1, -1, 1, 1)
    x <- add_baseline_data(x, data$baseline_data)
    expect_error(fit(x, refit = FALSE), "baseline|prepare_rtmle_data")
    x <- model_formula(prepare_rtmle_data(x), verbose = FALSE)
    x <- fit(x, refit = FALSE)
    fresh <- fit(continuous_analysis(data))

    expect_equal(x$estimate$Main_analysis, fresh$estimate$Main_analysis)
    expect_equal(x$continuous_nuisance$Always_A$data$L_0, data$baseline_data$L_0)
    expect_equal(x$continuous_nuisance$Always_A$data$propensity_0,
                 stats::plogis(data$baseline_data$L_0))
    expect_gt(abs(x$estimate$Main_analysis[Regime == "Always_A", Estimate] - before),
               0.05)
})

test_that("replacing a continuous regime invalidates its propensity fits", {
    data <- continuous_example()
    data$baseline_data$A_0 <- c(0, 0, 0, 1)
    data$timevarying_data$A <- data$baseline_data$A_0
    data$timevarying_data$event <- c("Y", "tauend", "tauend", "Y")
    x <- run_rtmle(continuous_analysis(data), estimator = "one-step")
    expect_equal(x$estimate$Main_analysis[Regime == "Always_A", IPW], 1,
                 tolerance = 1e-7)

    x <- regime(x, "Always_A", intervention = 0, treatment_variables = "A")
    x <- run_rtmle(x, estimator = "one-step", refit = FALSE)
    expect_equal(x$continuous_nuisance$Always_A$static_intervention, 0)
    expect_equal(x$continuous_nuisance$Always_A$data$propensity_0, rep(0.75, 4),
                 tolerance = 1e-7)
    expect_equal(x$estimate$Main_analysis[Regime == "Always_A", IPW], 1 / 3,
                 tolerance = 1e-7)
    expect_equal(x$estimate$Main_analysis[Regime == "Always_A", Estimate], 1 / 3,
                 tolerance = 1e-7)
})

test_that("changing estimator keeps estimates and influence curves consistent", {
    x <- run_rtmle(continuous_analysis(), estimator = "one-step")
    x <- run_rtmle(x, estimator = "tmle", refit = FALSE)
    expect_equal(as.character(unique(x$estimate$Main_analysis$Estimator)), "tmle")
    expect_equal(nrow(x$estimate$Main_analysis), 2L)
    expect_equal(mean(x$IC$risk$Always_A$time_horizon_1), 0, tolerance = 1e-8)
    x <- run_rtmle(x, estimator = "sequential_regression", refit = FALSE)
    expect_true(all(is.na(x$estimate$Main_analysis$Standard_error)))
    expect_true(all(is.na(summary(x)$Standard_error)))
})

test_that("continuous sequential regression uses censoring weights without fitting propensities", {
    data <- continuous_example()
    data$baseline_data$L_0 <- rep(0:1, 2)
    data$timevarying_data$L <- data$baseline_data$L_0
    data$timevarying_data$event[1L] <- "C"
    data$timevarying_data$time[1L] <- 1
    x <- continuous_analysis(data, regimes = "Always_A")

    outcome <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = rep(mean(data$rtmle_outcome), nrow(intervened_data)))
    }
    allow_propensity <- FALSE
    propensity_calls <- 0L
    propensity <- function(character_formula, data, intervened_data, ...) {
        if (!allow_propensity) stop("Sequential regression must not fit treatment propensities.")
        propensity_calls <<- propensity_calls + 1L
        list(predicted_values = rep(0.5, nrow(intervened_data)))
    }
    censoring <- list(name = "cox", fun = "learn_coxph", marginal_censoring_formula = ~ 1)
    fit <- function(x, estimator, ...) {
        run_rtmle(x, estimator = estimator, learner = outcome,
                  propensity_learner = propensity, censoring_learner = censoring, ...)
    }

    # At time 1 the Cox censoring hazard jumps by 1/4. The outcome at
    # time 2 contributes exp(1/4) to the pseudo-outcome regression.
    pooled <- fit(x, "sequential_regression")
    estimate <- pooled$estimate$Main_analysis
    expect_identical(as.character(estimate$Estimator), "sequential_regression")
    expect_equal(estimate$Estimate, exp(1/4) / 4)
    expect_equal(estimate$ICE_IPCW, estimate$Estimate)
    expect_identical(estimate$IPW, NA_real_)
    expect_identical(estimate$Standard_error, NA_real_)
    expect_equal(pooled$IC$risk$Always_A$time_horizon_1, rep(NA_real_, 4))
    expect_length(pooled$continuous_nuisance$Always_A$models$treatment, 0L)

    # Only the two baseline-treated people enter the stratified regression.
    stratified <- fit(pooled, "g-formula", refit = FALSE,
                      stratify_sequential_regression = TRUE)
    expect_identical(as.character(stratified$estimate$Main_analysis$Estimator),
                     "sequential_regression")
    expect_equal(stratified$estimate$Main_analysis$Estimate, exp(1/4) / 2)

    allow_propensity <- TRUE
    updated <- fit(pooled, "one-step", refit = FALSE)
    expect_identical(propensity_calls, 1L)
    expect_equal(updated$continuous_nuisance$Always_A$data$propensity_0, rep(0.5, 4))
})
