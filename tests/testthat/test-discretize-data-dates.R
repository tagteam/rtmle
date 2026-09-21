library(testthat)
library(data.table)
library(rtmle)


make_calendar_date_object <- function() {
    starts <- as.Date(c("2020-01-01", "2020-01-03", "2020-01-05"))
    x <- rtmle_init(
        time_grid = c(0, 10, 20),
        name_id = "id",
        name_outcome = "Y",
        name_competing = "D",
        name_censoring = "C",
        censored_label = "censored"
    )
    x <- add_long_data(
        x,
        outcome_data = data.table(
            id = 1:3,
            date = starts + c(12, 18, 15)
        ),
        competing_data = data.table(
            id = c(1, 3),
            date = starts[c(1, 3)] + c(20, 10)
        ),
        censored_data = data.table(
            id = 2,
            date = starts[2] + 5
        ),
        timevar_data = list(
            A = data.table(
                id = c(1, 1, 2, 3),
                date = starts[c(1, 1, 2, 3)] + c(0, 11, 0, 9),
                value = c(1, 0, 1, 1)
            ),
            exposure = data.table(
                id = 1:3,
                start_date = starts + 2,
                end_date = starts + 12
            )
        )
    )
    x <- add_baseline_data(
        x,
        data.table(id = 1:3, start = starts, age = c(50, 60, 70))
    )
    x
}


test_that("discretize_data converts calendar dates to numeric time since follow-up", {
    x <- make_calendar_date_object()
    original_starts <- x$data$baseline_data$start

    x <- discretize_data(
        x,
        start_followup_date = "start",
        A = "locf",
        exposure = "exposure_percent",
        verbose = FALSE
    )

    expect_true(is.numeric(x$long_data$outcome_data$date))
    expect_true(is.numeric(x$long_data$competing_data$date))
    expect_true(is.numeric(x$long_data$censored_data$date))
    expect_true(is.numeric(x$long_data$timevar_data$A$date))
    expect_true(is.numeric(x$long_data$timevar_data$exposure$start_date))
    expect_true(is.numeric(x$long_data$timevar_data$exposure$end_date))

    expect_equal(x$long_data$outcome_data$date, c(12, 18, 15))
    expect_equal(x$long_data$competing_data$date, c(20, 10))
    expect_equal(x$long_data$censored_data$date, 5)
    expect_equal(x$long_data$timevar_data$A$date, c(0, 11, 0, 9))
    expect_equal(x$long_data$timevar_data$exposure$start_date, c(2, 2, 2))
    expect_equal(x$long_data$timevar_data$exposure$end_date, c(12, 12, 12))
    expect_true(isTRUE(x$progress$substracted_start_followup_date))
    expect_s3_class(original_starts, "Date")
    expect_s3_class(x$data$baseline_data$start, "Date")

    expect_true(is.numeric(x$data$outcome_data$Y_0))
    expect_true(is.numeric(x$data$timevar_data$A$A_0))
    expect_true(is.numeric(x$data$timevar_data$exposure$exposure_1))

    # Calling discretize_data again must not subtract the subject-specific start
    # date a second time.
    outcome_dates <- x$long_data$outcome_data$date
    x <- discretize_data(
        x,
        start_followup_date = "start",
        A = "locf",
        exposure = "exposure_percent",
        verbose = FALSE
    )
    expect_equal(x$long_data$outcome_data$date, outcome_dates)
})


test_that("calendar-date histories require a baseline follow-up date", {
    x <- make_calendar_date_object()
    expect_error(
        discretize_data(x, A = "locf", verbose = FALSE),
        "Calendar-date.*require.*start_followup_date"
    )
})


test_that("adding long data resets the date-conversion state", {
    x <- make_calendar_date_object()
    x <- discretize_data(
        x,
        start_followup_date = "start",
        A = "locf",
        exposure = "exposure_percent",
        verbose = FALSE
    )
    x <- add_long_data(
        x,
        outcome_data = data.table::copy(x$long_data$outcome_data)
    )
    expect_false(isTRUE(x$progress$substracted_start_followup_date))
})


test_that("long_to_wide remains an obsolete alias", {
    x <- make_calendar_date_object()
    expect_warning(
        y <- long_to_wide(
            x,
            start_followup_date = "start",
            A = "locf",
            exposure = "exposure_percent",
            verbose = FALSE
        ),
        "deprecated"
    )
    expect_true(is.numeric(y$data$outcome_data$Y_0))
})


make_baseline_lookback_object <- function() {
    x <- rtmle_init(
        time_grid = c(0, 10, 20),
        name_id = "id",
        name_outcome = "Y",
        name_competing = "D",
        name_censoring = "C",
        censored_label = "censored"
    )
    x <- add_long_data(
        x,
        outcome_data = data.table(
            id = 1:2,
            date = c(15, 15)
        ),
        competing_data = data.table(
            id = 1,
            date = 18
        ),
        censored_data = data.table(
            id = 2,
            date = 12
        ),
        timevar_data = list(
            bleeding = data.table(
                id = c(1, 2),
                date = c(-5, 5)
            ),
            exposure = data.table(
                id = c(1, 2),
                start_date = c(-5, 5),
                end_date = c(0, 15)
            ),
            measurement = data.table(
                id = c(1, 2),
                date = c(-5, 5),
                value = c(10, 20)
            )
        )
    )
    add_baseline_data(x, data.table(id = 1:2))
}


test_that("baseline_lookback extends only the time-varying grid", {
    x_without_lookback <- discretize_data(
        make_baseline_lookback_object(),
        bleeding = "event_interval",
        exposure = "any_exposure",
        measurement = list(method = "measurement", fun_aggregate = "mean"),
        verbose = FALSE
    )
    x_with_lookback <- discretize_data(
        make_baseline_lookback_object(),
        bleeding = "event_interval",
        exposure = "any_exposure",
        measurement = list(method = "measurement", fun_aggregate = "mean"),
        baseline_lookback = 10,
        verbose = FALSE
    )

    # The pre-baseline bleeding and exposure are included in interval 0 only
    # when the time-varying grid is extended.
    expect_equal(
        x_without_lookback$data$timevar_data$bleeding$bleeding_0,
        c(0, 0)
    )
    expect_equal(
        x_with_lookback$data$timevar_data$bleeding$bleeding_0,
        c(1, 0)
    )
    expect_equal(
        x_without_lookback$data$timevar_data$exposure$exposure_0,
        c(0, 0)
    )
    expect_equal(
        x_with_lookback$data$timevar_data$exposure$exposure_0,
        c(1, 0)
    )
    expect_false(
        "measurement_0" %in%
            names(x_without_lookback$data$timevar_data$measurement)
    )
    expect_equal(
        x_with_lookback$data$timevar_data$measurement$measurement_0,
        c(10, NA_real_)
    )

    # Outcome histories are mapped with the unchanged follow-up grid.
    expect_identical(
        x_with_lookback$data$outcome_data,
        x_without_lookback$data$outcome_data
    )
})


test_that("baseline_lookback validates its value", {
    x <- make_baseline_lookback_object()
    expect_error(
        discretize_data(x, baseline_lookback = -1, verbose = FALSE),
        "single finite non-negative numeric"
    )
    expect_error(
        discretize_data(x, baseline_lookback = Inf, verbose = FALSE),
        "single finite non-negative numeric"
    )
    expect_error(
        discretize_data(x, baseline_lookback = "10", verbose = FALSE),
        "single finite non-negative numeric"
    )
})


make_baseline_exposure_start_object <- function() {
    x <- rtmle_init(
        time_grid = c(0, 10, 20),
        name_id = "id",
        name_outcome = "Y"
    )
    x <- add_long_data(
        x,
        outcome_data = data.table(
            id = 1:2,
            date = c(15, 15)
        ),
        timevar_data = list(
            treatment = data.table(
                id = 1:2,
                start_date = c(0, 5),
                end_date = c(10, 15)
            )
        )
    )
    add_baseline_data(x, data.table(id = 1:2))
}


test_that("baseline_exposure_start controls treatment exposure at time zero", {
    exposure_methods <- c(
        "exposure_time",
        "exposure_percent",
        "any_exposure",
        "has_exposure"
    )

    for (method in exposure_methods) {
        x_with_start_override <- discretize_data(
            make_baseline_exposure_start_object(),
            treatment = list(
                method = method,
                baseline_exposure_start = TRUE
            ),
            baseline_lookback = 10,
            verbose = FALSE
        )
        x_without_start_override <- discretize_data(
            make_baseline_exposure_start_object(),
            treatment = list(
                method = method,
                baseline_exposure_start = FALSE
            ),
            baseline_lookback = 10,
            verbose = FALSE
        )

        expect_equal(
            x_with_start_override$data$timevar_data$treatment$treatment_0,
            c(1, 0),
            info = method
        )
        expect_equal(
            x_without_start_override$data$timevar_data$treatment$treatment_0,
            c(0, 0),
            info = method
        )
    }
})


test_that("baseline_exposure_start must be logical", {
    expect_error(
        discretize_data(
            make_baseline_exposure_start_object(),
            treatment = list(
                method = "exposure_percent",
                baseline_exposure_start = 1
            ),
            baseline_lookback = 10,
            verbose = FALSE
        ),
        "single TRUE/FALSE"
    )
})
