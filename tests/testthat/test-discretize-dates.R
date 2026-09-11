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


test_that("discretize converts calendar dates to numeric time since follow-up", {
    x <- make_calendar_date_object()
    original_starts <- x$data$baseline_data$start

    x <- discretize(
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

    # Calling discretize again must not subtract the subject-specific start
    # date a second time.
    outcome_dates <- x$long_data$outcome_data$date
    x <- discretize(
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
        discretize(x, A = "locf", verbose = FALSE),
        "Calendar-date.*require.*start_followup_date"
    )
})


test_that("adding long data resets the date-conversion state", {
    x <- make_calendar_date_object()
    x <- discretize(
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
