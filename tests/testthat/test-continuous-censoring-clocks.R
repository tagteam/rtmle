# Two horizons retain different last event indices. At the shorter horizon,
# some subjects' next events are pooled into their eventual terminal event.
continuous_clock_fixture <- function() {
    baseline <- data.frame(id = 1:6, A_0 = c(0, 1, 0, 1, 0, 1),
                           L_0 = 1:6, age = 41:46)
    one_person <- function(id, times, events, treatment) {
        data.frame(id = id, time = times, event = events, A = treatment,
                   L = id + seq_along(times) / 10)
    }
    events <- do.call(rbind, list(
        one_person(1, c(20, 45, 95), c("A", "L", "C"), 1),
        one_person(2, c(80, 90, 120), c("A", "L", "Y"), 0),
        one_person(3, c(30, 70), c("A", "Y"), 1),
        one_person(4, c(15, 50, 130), c("A", "L", "Y"), 0),
        one_person(5, c(25, 60), c("A", "C"), 1),
        one_person(6, c(35, 65, 75, 201), c("A", "L", "L", "tauend"), 0)
    ))
    prepare_continuous_data(list(baseline_data = baseline, timevarying_data = events),
                            time_horizons = c(100, 200),
                            time_covariates = c("A", "L"),
                            baseline_covariates = c("A_0", "L_0", "age"),
                            min_events = 0)
}

continuous_clock_spy <- function() {
    calls <- list()
    learner <- function(character_formula, data, intervened_data, time_variable,
                        time_ref = NULL, reuse_fit = NULL, ...) {
        calls[[length(calls) + 1L]] <<- list(
            time_variable = time_variable, time_ref = time_ref,
            reused = !is.null(reuse_fit), id = intervened_data$id,
            endpoint = intervened_data[[time_variable]],
            start = if (is.null(time_ref)) NULL else intervened_data[[time_ref]])
        list(predicted_values = rep(0.8, nrow(intervened_data)), fit = list())
    }
    list(learner = learner, calls = function() calls)
}

continuous_clock_call <- function(calls, variable, reused) {
    selected <- Filter(function(call) identical(call$time_variable, variable) &&
                           identical(call$reused, reused), calls)
    if (length(selected) != 1L) stop("Expected one censoring call for ", variable)
    selected[[1L]]
}

continuous_clock_value <- function(call, id, field) {
    call[[field]][match(id, call$id)]
}

# Propensity estimates are immaterial here; the recording censoring learner
# exposes the actual times passed to each censoring fit and prediction.
continuous_clock_treatment <- function(character_formula, data, intervened_data, ...) {
    list(predicted_values = rep(0.5, nrow(intervened_data)))
}

test_that("marginal censoring uses absolute times and conditions on the previous event", {
    prepared <- continuous_clock_fixture()
    expect_equal(prepared$info$last_event, c(2, 3))
    spy <- continuous_clock_spy()
    continuous_propensity(prepared, learner = continuous_clock_treatment,
                          censoring_learner = list(name = "clock", fun = spy$learner,
                                                   marginal_censoring = TRUE))
    calls <- spy$calls()
    initial <- continuous_clock_call(calls, "time", FALSE)
    expect_null(initial$time_ref)
    expect_equal(continuous_clock_value(initial, 1, "endpoint"), 95)
    expect_equal(continuous_clock_value(initial, 2, "endpoint"), 120)

    first <- continuous_clock_call(calls, "time_1", TRUE)
    expect_equal(continuous_clock_value(first, 1, "endpoint"), 20)
    expect_equal(continuous_clock_value(first, 1, "start"), 0)
    ordinary <- continuous_clock_call(calls, "time_2", TRUE)
    pooled <- continuous_clock_call(calls, "time_pooled_2", TRUE)
    expect_identical(ordinary$time_ref, "time_1")
    expect_identical(pooled$time_ref, "time_1")
    expect_equal(continuous_clock_value(ordinary, 1, "endpoint"), 45)
    expect_equal(continuous_clock_value(pooled, 1, "endpoint"), 95)
    expect_equal(continuous_clock_value(pooled, 1, "start"), 20)
    expect_equal(continuous_clock_value(pooled, 2, "endpoint"), 120)
    expect_equal(continuous_clock_value(pooled, 2, "start"), 80)

    # The maximum-index pooled endpoint can differ from the ordinary endpoint.
    ordinary_max <- continuous_clock_call(calls, "time_3", TRUE)
    pooled_max <- continuous_clock_call(calls, "time_pooled_3", TRUE)
    expect_equal(continuous_clock_value(ordinary_max, 6, "endpoint"), 75)
    expect_equal(continuous_clock_value(pooled_max, 6, "endpoint"), 201)
    expect_equal(continuous_clock_value(pooled_max, 6, "start"), 65)
})

test_that("conditional censoring uses elapsed times for ordinary and pooled events", {
    prepared <- continuous_clock_fixture()
    spy <- continuous_clock_spy()
    continuous_propensity(prepared, learner = continuous_clock_treatment,
                          censoring_learner = list(name = "clock", fun = spy$learner,
                                                   marginal_censoring = FALSE))
    calls <- spy$calls()
    expect_true(all(vapply(calls, function(call) is.null(call$time_ref) && !call$reused,
                           logical(1))))
    first <- continuous_clock_call(calls, "time_1", FALSE)
    expect_equal(continuous_clock_value(first, 1, "endpoint"), 20)
    ordinary <- continuous_clock_call(calls, "time_2", FALSE)
    pooled <- continuous_clock_call(calls, "time_pooled_2", FALSE)
    expect_equal(continuous_clock_value(ordinary, 1, "endpoint"), 25)
    expect_equal(continuous_clock_value(pooled, 1, "endpoint"), 75)
    expect_equal(continuous_clock_value(ordinary, 2, "endpoint"), 10)
    expect_equal(continuous_clock_value(pooled, 2, "endpoint"), 40)

    ordinary_max <- continuous_clock_call(calls, "time_3", FALSE)
    pooled_max <- continuous_clock_call(calls, "time_pooled_3", FALSE)
    expect_equal(continuous_clock_value(ordinary_max, 6, "endpoint"), 10)
    expect_equal(continuous_clock_value(pooled_max, 6, "endpoint"), 136)
})

test_that("Cox and Weibull learners use conditional survival over absolute times", {
    training <- data.frame(time = c(5, 20, 45, 75, 95, 130),
                           status = c(1, 1, 1, 1, 1, 0))
    prediction <- data.frame(end = 95, start = 20)
    cox <- learn_coxph("survival::Surv(time, status) ~ 1", training, prediction,
                       time_variable = "end", time_ref = "start")
    hazard <- cox$baseline_hazard
    H_before_end <- tail(c(0, hazard$hazard[hazard$time < 95]), 1)
    H_at_start <- tail(c(0, hazard$hazard[hazard$time <= 20]), 1)
    expect_equal(cox$predicted_values, exp(-(H_before_end - H_at_start)))

    weibull <- learn_weibull("survival::Surv(time, status) ~ 1", training, prediction,
                             time_variable = "end", time_ref = "start")
    location <- unname(stats::predict(weibull$fit, newdata = prediction, type = "lp"))
    shape <- 1 / weibull$fit$scale
    H <- function(t) exp(shape * (log(t) - location))
    expect_equal(weibull$predicted_values, exp(-(H(95) - H(20))))
})
