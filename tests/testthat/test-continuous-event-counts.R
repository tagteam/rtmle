event_count_analysis <- function(baseline_counts = NULL) {
    baseline <- data.frame(id = 1:2, age = c(50, 60))
    if (!is.null(baseline_counts)) baseline$bleeding_0 <- baseline_counts
    x <- rtmle_init(c(0, 5), name_id = "id", name_time = "date",
                    name_outcome = "Y", continuous_time = TRUE)
    x <- add_baseline_data(x, baseline)
    x <- add_continuous_data(
        x, outcome_data = data.frame(id = 1, date = 5),
        censored_data = data.frame(id = 2, date = 5, event = "tauend"),
        timevar_data = list(
            A = data.frame(id = rep(1:2, each = 3), date = rep(c(0, 2, 4), 2),
                           value = c(0, 0, 1, 1, 1, 0)),
            bleeding = data.frame(id = c(1, 1, 1, 2, 2, 2),
                                  date = c(1, 2, 3, 1, 1, 3))),
        treatment_events = "A", covariate_events = "bleeding")
    prepare_rtmle_data(x, min_events = 0)
}

test_that("prepared histories carry cumulative event counts through later events", {
    x <- event_count_analysis()
    expect_equal(x$prepared_data$bleeding_0, c(0, 0))
    expect_equal(x$prepared_data$bleeding_1, c(1, 1))
    expect_equal(x$prepared_data$bleeding_2, c(2, 1))
    expect_equal(x$prepared_data$bleeding_3, c(3, 2))
    expect_equal(x$prepared_data$bleeding_4, c(3, 2))
    expect_equal(x$prepared_data$bleeding_5, c(3, 2))
    # The simultaneous event and treatment opportunity is one A-marked event.
    expect_equal(as.character(x$prepared_data$event_2), c("A", "A"))
})

test_that("factor baseline counts become numeric and increment from their values", {
    x <- event_count_analysis(factor(c(2, 4)))
    expect_type(x$prepared_data$bleeding_0, "double")
    expect_equal(x$prepared_data$bleeding_0, c(2, 4))
    expect_equal(x$prepared_data$bleeding_1, c(3, 5))
    expect_equal(x$prepared_data$bleeding_3, c(5, 6))
})
