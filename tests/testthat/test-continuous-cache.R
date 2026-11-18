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
    x <- run_rtmle(x, estimator = "g-formula", refit = FALSE)
    expect_true(all(is.na(x$estimate$Main_analysis$Standard_error)))
    expect_true(all(is.na(summary(x)$Standard_error)))
})
