# Administrative follow-up is inferred after considering every terminal stream.
tauend_fixture <- function(horizons = c(2, 5)) {
    baseline <- data.frame(person = 1:6, A = 0)
    x <- rtmle_init(c(0, horizons), name_id = "person", name_time = "date",
                    name_outcome = "stroke", name_competing = "death",
                    name_censoring = "dropout", continuous_time = TRUE)
    add_baseline_data(x, baseline)
}

tauend_add_events <- function(x, censored = NULL, outcome = NULL, competing = NULL) {
    args <- list(x = x, treatment_events = "A", covariate_events = NULL)
    if (!is.null(censored)) args$censored_data <- censored
    if (!is.null(outcome)) args$outcome_data <- outcome
    if (!is.null(competing)) args$competing_data <- competing
    do.call(add_continuous_data, args)
}

tauend_first <- function(x) {
    out <- as.data.frame(x$prepared_data)
    out <- out[order(out$person), ]
    data.frame(person = out$person, time = out$time_1,
               event = as.character(out$event_1))
}

tauend_censor <- function(x) {
    out <- as.data.frame(x$long_data$censored_data)
    out <- out[order(out$person), c("person", "date", "event"), drop = FALSE]
    rownames(out) <- NULL
    out$event <- as.character(out$event)
    out
}

test_that("two-column censoring records are C and missing terminal records end one unit beyond the horizon", {
    x <- tauend_fixture()
    censor <- data.frame(person = c(1, 2), date = c(3, 6))
    outcome <- data.frame(person = 3, date = 4)
    competing <- data.frame(person = 4, date = 4.5)
    before_censor <- data.table::copy(censor)
    before_outcome <- data.table::copy(outcome)
    before_competing <- data.table::copy(competing)
    x <- tauend_add_events(x, censor, outcome, competing)
    result <- prepare_rtmle_data(x)

    expect_equal(tauend_first(result), data.frame(
        person = 1:6, time = c(3, 6, 4, 4.5, 6, 6),
        event = c("C", "C", "Y", "D", "tauend", "tauend")))
    # Censoring beyond the horizon is still a C event; it is simply not
    # censoring before that horizon.
    expect_equal(result$continuous_prepared$info$is_censored, c(FALSE, TRUE))
    expect_equal(tauend_censor(result), data.frame(
        person = c(1, 2, 5, 6), date = c(3, 6, 6, 6),
        event = c("C", "C", "tauend", "tauend")))
    expect_identical(censor, before_censor)
    expect_identical(outcome, before_outcome)
    expect_identical(competing, before_competing)
})

test_that("tauend is inferred without a censoring table and depends on the initialized horizon", {
    for (horizon in c(5, 8)) {
        x <- tauend_fixture(horizons = c(2, horizon))
        x <- tauend_add_events(x,
            outcome = data.frame(person = 1, date = 1),
            competing = data.frame(person = 2, date = 1.5))
        result <- prepare_rtmle_data(x)
        expected <- data.frame(person = 1:6,
                               time = c(1, 1.5, rep(horizon + 1, 4)),
                               event = c("Y", "D", rep("tauend", 4)))
        expect_equal(tauend_first(result), expected)
        expected_censor <- data.frame(person = 3:6, date = rep(horizon + 1, 4),
                                      event = rep("tauend", 4))
        expect_equal(tauend_censor(result), expected_censor)
        repeated <- prepare_rtmle_data(result)
        expect_equal(tauend_first(repeated), expected)
        expect_equal(tauend_censor(repeated), expected_censor)
    }
})

test_that("explicit administrative events are preserved", {
    x <- tauend_fixture()
    censor <- data.frame(person = c(1, 2), date = c(6, 7), event = c("C", "tauend"))
    result <- prepare_rtmle_data(tauend_add_events(x, censor))
    expect_equal(tauend_first(result), data.frame(
        person = 1:6, time = c(6, 7, 6, 6, 6, 6),
        event = c("C", "tauend", rep("tauend", 4))))
    expect_false(any(result$continuous_prepared$info$is_censored))
    expect_equal(tauend_censor(result), data.frame(
        person = 1:6, date = c(6, 7, 6, 6, 6, 6),
        event = c("C", "tauend", rep("tauend", 4))))
})

test_that("calendar terminal dates are rejected", {
    x <- tauend_add_events(tauend_fixture(),
                           censored = data.frame(person = 1, date = as.Date("2020-01-03")))
    expect_error(prepare_rtmle_data(x), "numeric.*since baseline")
})

test_that("removed preparation options are rejected", {
    x <- tauend_add_events(tauend_fixture())
    expect_error(prepare_rtmle_data(x, start_followup_date = 0), "start_followup_date")
    expect_error(prepare_rtmle_data(x, baseline_lookback = 0), "baseline_lookback")
})

test_that("terminal records after the horizon remain their recorded event type", {
    x <- tauend_add_events(tauend_fixture(),
        censored = data.frame(person = 1, date = 5),
        outcome = data.frame(person = 2, date = 7),
        competing = data.frame(person = 3, date = 8))
    result <- prepare_rtmle_data(x)
    expect_equal(tauend_first(result), data.frame(
        person = 1:6, time = c(5, 7, 8, 6, 6, 6),
        event = c("C", "Y", "D", rep("tauend", 3))))
    expect_equal(result$continuous_prepared$info$is_censored, c(FALSE, TRUE))
})

test_that("all subjects without terminal data receive inferred administrative endpoints", {
    result <- prepare_rtmle_data(tauend_add_events(tauend_fixture()))
    expect_equal(tauend_first(result), data.frame(
        person = 1:6, time = rep(6, 6), event = rep("tauend", 6)))
    expect_false(any(result$continuous_prepared$info$is_censored))
})

test_that("existing outcome and competing events retain earliest-event precedence", {
    x <- tauend_add_events(tauend_fixture(),
        censored = data.frame(person = 1:3, date = c(3, 2, 5)),
        outcome = data.frame(person = c(1, 3), date = c(2, 4)),
        competing = data.frame(person = 2, date = 4))
    result <- prepare_rtmle_data(x)
    expect_equal(tauend_first(result), data.frame(
        person = 1:6, time = c(2, 2, 4, 6, 6, 6),
        event = c("Y", "C", "Y", rep("tauend", 3))))
})

test_that("tied outcome and censoring events are rejected", {
    x <- tauend_add_events(tauend_fixture(),
        censored = data.frame(person = 3, date = 4),
        outcome = data.frame(person = 3, date = 4))
    expect_error(prepare_rtmle_data(x), "[Tt]ie|[Tt]ied|same time")
})

test_that("inferred censoring rows are recalculated after changing the observed data", {
    x <- prepare_rtmle_data(tauend_add_events(tauend_fixture()))
    expect_equal(tauend_censor(x), data.frame(person = 1:6, date = rep(6, 6),
                                                event = rep("tauend", 6)))

    # Adding a real outcome must remove that subject's inferred record.
    x <- add_continuous_data(x, outcome_data = data.frame(person = 1, date = 3))
    x <- prepare_rtmle_data(x)
    expect_equal(tauend_censor(x), data.frame(person = 2:6, date = rep(6, 5),
                                               event = rep("tauend", 5)))
    expect_equal(tauend_first(x)$event, c("Y", rep("tauend", 5)))

    # Adding censoring replaces inferred records, and a later preparation
    # restores only those still needed for subjects without a terminal event.
    x <- add_continuous_data(x, censored_data = data.frame(person = 2, date = 4))
    x <- prepare_rtmle_data(x)
    expect_equal(tauend_censor(x), data.frame(person = 2:6, date = c(4, rep(6, 4)),
                                               event = c("C", rep("tauend", 4))))
    expect_equal(tauend_first(x)$event, c("Y", "C", rep("tauend", 4)))

    # A changed baseline cohort must not retain an inferred endpoint for a
    # removed subject and must infer one for a newly added subject.
    x <- add_baseline_data(x, data.frame(person = c(1:5, 7), A = 0))
    x <- prepare_rtmle_data(x)
    expect_equal(tauend_censor(x), data.frame(person = c(2:5, 7), date = c(4, rep(6, 4)),
                                               event = c("C", rep("tauend", 4))))
    expect_equal(tauend_first(x)$person, c(1:5, 7))
    expect_equal(tauend_first(x)$event, c("Y", "C", rep("tauend", 4)))
})

test_that("a replacement tauend matching an inferred row remains explicit", {
    for (outcome_together in c(FALSE, TRUE)) {
        x <- prepare_rtmle_data(tauend_add_events(tauend_fixture()))
        explicit <- data.frame(person = 1, date = 6, event = "tauend")
        outcome <- data.frame(person = 1, date = 3)
        x <- if (outcome_together) {
            add_continuous_data(x, censored_data = explicit, outcome_data = outcome)
        } else {
            add_continuous_data(x, censored_data = explicit)
        }
        x <- prepare_rtmle_data(x)
        expect_false(1 %in% x$continuous_inferred_censoring$person)

        if (!outcome_together) {
            x <- add_continuous_data(x, outcome_data = outcome)
            x <- prepare_rtmle_data(x)
        }
        expect_equal(tauend_censor(x), data.frame(person = 1:6, date = rep(6, 6),
                                                  event = rep("tauend", 6)))
        expect_equal(tauend_first(x)$event[1], "Y")
    }
})

test_that("fractional horizons retain their scale", {
    x <- tauend_fixture(horizons = c(1, 2.5))
    x <- tauend_add_events(x, data.frame(person = 1, date = 1.5))
    result <- prepare_rtmle_data(x)
    expect_equal(tauend_first(result), data.frame(person = 1:6,
        time = c(1.5, rep(3.5, 5)), event = c("C", rep("tauend", 5))))
    expect_equal(tauend_censor(result), data.frame(person = 1:6,
        date = c(1.5, rep(3.5, 5)), event = c("C", rep("tauend", 5))))
})

