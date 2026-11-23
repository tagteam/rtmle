test_that("continuous analyses reuse the rtmle pipeline and keep the baseline TMLE tilt", {
    x <- continuous_analysis()
    constant <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = rep(0.2, nrow(intervened_data)))
    }
    x <- run_rtmle(x, learner = list(name = "constant", fun = constant),
                   propensity_learner = "glm", estimator = "tmle")
    e <- x$estimate$Main_analysis
    expect_equal(e[Regime == "Always_A", Estimate], 0.5, tolerance = 1e-7)
    expect_equal(e[Regime == "Never_A", Estimate], 0, tolerance = 1e-7)
    expect_equal(mean(x$IC$risk$Always_A$time_horizon_1), 0, tolerance = 1e-7)
    expect_s3_class(summary(x), "data.table")
    expect_true(nrow(summary(x)) > 0)
    expect_equal(summary(x)[Target_parameter == "Risk_difference", Estimate], -0.5)
})

test_that("risk plots use original continuous time units", {
    d <- continuous_example()
    d$timevarying_data$time <- 720
    x <- run_rtmle(continuous_analysis(d, horizons = c(100, 720)), estimator = "one-step")
    p <- ggplot2::autoplot(x)
    expect_equal(p$scales$get_scales("x")$breaks, c(0, 1, 2))
    expect_equal(p$scales$get_scales("x")$labels, c("0", "100", "720"))
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("both static regimes have the correct adherence weights", {
    x <- run_rtmle(continuous_analysis(), estimator = "one-step")
    e <- x$estimate$Main_analysis
    expect_equal(e[Regime == "Always_A", IPW], 0.5)
    expect_equal(e[Regime == "Never_A", IPW], 0)
    expect_equal(e[Regime == "Always_A", Estimate], 0.5, tolerance = 1e-7)
    expect_equal(e[Regime == "Never_A", Estimate], 0, tolerance = 1e-7)
    expect_equal(length(x$IC$risk$Always_A$time_horizon_1), 4L)
})

test_that("preparation copies and orders data and handles an empty outcome history", {
    d <- continuous_example()
    before <- data.table::copy(d)
    p <- prepare_continuous_data(d, c(1, 2), c("A", "L"), c("A_0", "L_0"))
    expect_identical(d, before)
    expect_equal(p$info$last_event, c(1, 1))
    x <- run_rtmle(continuous_analysis(d, c(1, 2)), estimator = "one-step", time_horizon = c(1, 2))
    expect_equal(x$estimate$Main_analysis[Time_horizon == 1, Estimate], c(0, 0))
    expect_equal(x$estimate$Main_analysis[Time_horizon == 1, Standard_error], c(0, 0))
    d$timevarying_data$event <- "tauend"
    x <- run_rtmle(continuous_analysis(d), estimator = "tmle")
    expect_equal(x$estimate$Main_analysis$Estimate, c(0, 0))
})

test_that("invalid continuous input is rejected before model fitting", {
    d <- continuous_example()
    prep <- function(d, ...) prepare_continuous_data(d, 2, c("A", "L"), c("A_0", "L_0"), ...)
    bad <- d
    bad$baseline_data$id[2] <- 1
    expect_error(prep(bad), "unique")
    bad <- d
    bad$timevarying_data$id[1] <- 10
    expect_error(prep(bad), "same subjects")
    bad <- d
    bad$timevarying_data$A[1] <- 2
    expect_error(prep(bad), "binary")
    bad <- d
    bad$timevarying_data$L[1] <- NA
    expect_error(prep(bad), "missing|NA")
    bad <- d
    bad$timevarying_data <- rbind(bad$timevarying_data, bad$timevarying_data[1, ])
    expect_error(prep(bad), "ties")
    bad <- d
    bad$timevarying_data$event[1] <- "invalid"
    expect_error(prep(bad), "event")
    bad <- d
    bad$timevarying_data$time[1] <- Inf
    expect_error(prep(bad), "finite")
    bad <- d
    bad$timevarying_data$event[1] <- "A"
    expect_error(prep(bad), "terminal")
    expect_error(prepare_continuous_data(d, NA_real_, c("A", "L"), c("A_0", "L_0")), "finite")
    expect_error(prep(d, min_events = -1), "min_events")
    expect_error(run_rtmle(continuous_analysis(), conservative = FALSE, estimator = "tmle"), "one-step")
    expect_error(run_rtmle(continuous_analysis(), estimator = "unknown"), "estimator")
    expect_error(run_rtmle(continuous_analysis(), time_horizon = 3), "time_horizon")
    expect_error(regime(continuous_analysis(), "late", data.frame(time = 999, A = 1)), "scheduled")
    bad <- d
    bad$baseline_data$time_0 <- 1:4
    expect_error(prepare_continuous_data(bad, 2, c("A", "L"), c("A_0", "time_0")), "collide")
})

test_that("baseline covariate names need not mirror later measurements", {
    d <- continuous_example()
    names(d$baseline_data)[3] <- "initial_L"
    p <- prepare_continuous_data(d, 2, c("A", "L"), c("A_0", "initial_L"))
    expect_equal(nrow(p$wide_data), 4L)
})

test_that("outcome learner warnings remain available as diagnostics", {
    learner <- function(character_formula, data, intervened_data, ...) {
        warning("Outcome fit did not converge")
        list(predicted_values = rep(0.2, nrow(intervened_data)))
    }
    x <- run_rtmle(continuous_analysis(), learner = learner,
                   propensity_learner = "glm", estimator = "one-step")
    expect_true(any(grepl("Outcome fit did not converge", x$diagnostics$learner_warnings$warning)))
})

test_that("custom identifiers, times and event marks map without collisions", {
    d <- continuous_example()
    names(d$baseline_data)[1] <- "person"
    names(d$timevarying_data)[1:2] <- c("person", "date")
    d$timevarying_data$event[d$timevarying_data$event == "Y"] <- "C"
    x <- rtmle_init(c(0, 2), "person", name_time = "date", name_outcome = "C",
                    name_censoring = "Y", continuous_time = TRUE)
    x <- add_baseline_data(x, d$baseline_data)
    events <- d$timevarying_data
    x <- add_continuous_data(x,
        outcome_data = events[events$event == "C", c("person", "date")],
        censored_data = events[events$event == "tauend", c("person", "date", "event")],
        timevar_data = list(A = data.frame(person = integer(), date = numeric(), value = numeric()),
                             L = data.frame(person = integer(), date = numeric(), value = numeric())),
        treatment_events = "A", covariate_events = "L")
    expect_equal(names(d$timevarying_data)[1:2], c("person", "date"))
    x <- prepare_rtmle_data(x)
    expect_equal(sum(x$continuous_prepared$wide_data$event_1 == "Y"), 1L)
    expect_equal(sum(x$continuous_prepared$wide_data$event_1 == "C"), 0L)
    expect_true("person" %in% names(x$prepared_data))
})

test_that("continuous analyses use the same horizon indices as discrete analyses", {
    d <- continuous_example()
    d$timevarying_data$time <- 720
    x <- continuous_analysis(d, horizons = c(100, 720), regimes = "Always_A")
    first <- run_rtmle(x, estimator = "one-step", time_horizon = 1)
    expect_equal(first$estimate$Main_analysis$Time_horizon, 1)
    for (column in c("Target", "Regime", "Target_parameter", "Estimator"))
        expect_s3_class(first$estimate$Main_analysis[[column]], "factor")
    expect_true(all(is.na(first$estimate$Main_analysis$P_value)))
    expect_named(first$IC$risk$Always_A, "time_horizon_1")
    final <- run_rtmle(x, estimator = "one-step")
    expect_equal(final$estimate$Main_analysis$Time_horizon, 2)
    expect_named(final$IC$risk$Always_A, "time_horizon_2")
    expect_error(run_rtmle(x, estimator = "one-step", time_horizon = 720), "time_horizon")
})

test_that("continuous follow-up counts completed intervals at grid boundaries", {
    d <- continuous_example()
    d$timevarying_data$time <- c(30, 100, 101, 720)
    d$timevarying_data$event <- c("C", "Y", "C", "Y")
    x <- continuous_analysis(d, horizons = c(100, 720), regimes = "Always_A")
    expect_equal(x$followup$last_interval, c(0, 0, 1, 1))
})
