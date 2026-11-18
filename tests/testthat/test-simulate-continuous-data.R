test_that("continuous simulation returns ordered event and baseline data", {
    set.seed(34)
    simulated <- simulate_continuous_time_data(100)
    baseline <- simulated$baseline_data
    events <- simulated$timevarying_data
    expect_named(simulated, c("baseline_data", "timevarying_data"))
    expect_named(baseline, c("id", "age", "A_0", "L_0"))
    expect_named(events, c("id", "time", "event", "A", "L"))
    expect_equal(baseline$id, seq_len(100))
    expect_true(all(baseline$age >= 40 & baseline$age <= 90))
    expect_true(all(events$time > 0 & events$time <= 900))
    expect_true(all(events$event %in% c("A", "L", "C", "Y", "D", "tauend")))
    expect_true(all(events[, all(diff(time) > 0), by = id]$V1))
    terminal <- events[!event %in% c("A", "L")]
    expect_equal(terminal$id, seq_len(100))
})

test_that("continuous simulation supports no censoring and static interventions", {
    set.seed(11)
    simulated <- simulate_continuous_time_data(100, static_intervention = 0,
                                               no_competing_events = TRUE)
    expect_true(all(simulated$baseline_data$A_0 == 0))
    expect_true(all(simulated$timevarying_data$A == 0))
    expect_false(any(simulated$timevarying_data$event %in% c("C", "D")))

    simulated <- simulate_continuous_time_data(100, uncensored = TRUE,
                                               static_intervention_baseline = NULL)
    expect_false(any(simulated$timevarying_data$event == "C"))
    expect_setequal(simulated$baseline_data$A_0, c(0, 1))
})

test_that("continuous simulation handles zero rates and one event", {
    rates <- list(A = 0, L = 0, C = 0, Y = 0, D = 0)
    simulated <- simulate_continuous_time_data(5, K = 1,
                                               baseline_rate_list = rates,
                                               aggregate_parameters = TRUE,
                                               max_fup = 10)
    expect_equal(simulated$timevarying_data$time, rep(10, 5))
    expect_true(all(simulated$timevarying_data$event == "tauend"))
    expect_equal(nrow(simulated$baseline_data), 5)
})

test_that("continuous simulation validates sample size and event limits", {
    for (n in list(0, -1, 1.5, NA_real_, Inf, numeric(), c(1, 2), "5")) {
        expect_error(simulate_continuous_time_data(n), "n must")
    }
    for (K in list(0, -1, 1.5, NA_real_, Inf)) {
        expect_error(simulate_continuous_time_data(5, K = K), "K must")
    }
    expect_error(simulate_continuous_time_data(5, limit_event_A = -1),
                 "limit_event_A")
    expect_error(simulate_continuous_time_data(5, limit_event_L = 0.5),
                 "limit_event_L")
})

test_that("continuous simulation validates rates and follow-up settings", {
    expect_error(simulate_continuous_time_data(5, baseline_rate_list = list(A = 1)),
                 "baseline_rate_list")
    rates <- list(A = 0, L = 0, C = 0, Y = -1, D = 0)
    expect_error(simulate_continuous_time_data(5, baseline_rate_list = rates),
                 "baseline_rate_list")
    rates$Y <- NA_real_
    expect_error(simulate_continuous_time_data(5, baseline_rate_list = rates),
                 "baseline_rate_list")
    expect_error(simulate_continuous_time_data(5, max_fup = 0), "max_fup")
    expect_error(simulate_continuous_time_data(5, visitation_interval = -1),
                 "visitation_interval")
    expect_error(simulate_continuous_time_data(5, visitation_sd = -1),
                 "visitation_sd")
})

test_that("continuous simulation validates interventions, flags, and effects", {
    expect_error(simulate_continuous_time_data(5, static_intervention = 2),
                 "static_intervention")
    expect_error(simulate_continuous_time_data(5, static_intervention_baseline = 2),
                 "static_intervention_baseline")
    expect_error(simulate_continuous_time_data(5, static_intervention_stop = 1),
                 "static_intervention_stop.*not supported")
    expect_error(simulate_continuous_time_data(5, uncensored = NA), "uncensored")
    expect_error(simulate_continuous_time_data(5, no_competing_events = 1),
                 "no_competing_events")
    expect_error(simulate_continuous_time_data(5, effects = list()), "effects")
    expect_error(simulate_continuous_time_data(5, K = 4), "effects.*beta_")
    expect_no_error(simulate_continuous_time_data(5, K = 4,
                                                 aggregate_parameters = TRUE))
})

test_that("continuous simulation keeps scheduled visits forward in time", {
    set.seed(21)
    simulated <- simulate_continuous_time_data(100, visitation_interval = 1,
                                               visitation_sd = 100,
                                               uncensored = TRUE)
    expect_true(all(simulated$timevarying_data$time > 0))
    expect_true(all(simulated$timevarying_data[, all(diff(time) > 0), by = id]$V1))
})
