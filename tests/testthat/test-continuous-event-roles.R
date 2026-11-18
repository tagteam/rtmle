continuous_role_analysis <- function() {
    x <- rtmle_init(c(0, 3), "id", name_time = "date", name_outcome = "Y",
                    continuous_time = TRUE)
    x <- add_baseline_data(x, data.frame(id = 1:2, A = c(0, 1), L = 0))
    treatment <- data.frame(id = rep(1:2, each = 2), date = rep(c(1, 2), 2), value = c(0, 1, 1, 0))
    x <- add_continuous_data(x,
        outcome_data = data.frame(id = 1, date = 3),
        censored_data = data.frame(id = 2, date = 3, event = "tauend"),
        timevar_data = list(A = treatment, L = data.frame(id = 1:2, date = 1, value = c(1, 2)),
                            visits = data.frame(id = 1:2, date = 2)),
        treatment_events = "visits",
        covariate_events = "L")
    x
}

test_that("measurement of treatment does not itself create an opportunity", {
    x <- prepare_rtmle_data(continuous_role_analysis(), min_events = 0)
    expect_equal(as.character(x$prepared_data$event_1), c("L", "L"))
    expect_equal(as.character(x$prepared_data$event_2), c("A", "A"))
    expect_equal(x$prepared_data$A_1, c(0, 1))
    expect_equal(x$prepared_data$A_2, c(1, 0))
    expect_equal(x$prepared_data$L_1, c(1, 2))
    bad <- continuous_role_analysis()
    bad$long_data$timevar_data$A[1, value := 1]
    expect_error(prepare_rtmle_data(bad), "treatment")
})

test_that("covariate exposure intervals supply states at declared events", {
    x <- continuous_role_analysis()
    x <- add_continuous_data(x, B = data.frame(id = 1:2, start_date = 0.5, end_date = 1.5))
    before <- data.table::copy(x$long_data)
    prepared <- prepare_rtmle_data(x, min_events = 0)
    expect_equal(prepared$prepared_data$time_1, c(1, 1))
    expect_equal(prepared$prepared_data$time_2, c(2, 2))
    expect_equal(prepared$prepared_data$B_1, c(1, 1))
    expect_equal(prepared$prepared_data$B_2, c(0, 0))
    expect_equal(prepared$long_data, before)
    expect_equal(prepared$continuous_event_roles$treatment, "visits")
})

test_that("calendar event-role times are rejected", {
    x <- continuous_role_analysis()
    x$continuous_event_roles$covariate <- data.frame(id = 1, date = as.Date("2026-01-02"))
    expect_error(prepare_rtmle_data(x, min_events = 0), "numeric.*since baseline")
})

