continuous_treatment_input <- function() {
    x <- rtmle_init(c(0, 3), name_id = "id", name_time = "date",
                    name_outcome = "Y", name_censoring = "C",
                    continuous_time = TRUE)
    add_baseline_data(x, data.frame(id = 1:2, age = c(50, 60)))
}

continuous_treatment_points <- function() {
    data.frame(id = rep(1:2, each = 3), date = rep(0:2, 2),
               value = c(0, 0, 1, 1, 1, 0))
}

test_that("treatment event declarations are nonempty unique stream names", {
    x <- continuous_treatment_input()
    points <- continuous_treatment_points()
    invalid <- list(NULL, NA_character_, character(), "", c("A", "A"),
                    data.frame(id = 1, date = 1), 1)
    for (value in invalid) {
        expect_error(add_continuous_data(x, timevar_data = list(A = points),
                                         treatment_events = value,
                                         covariate_events = NULL),
                     "treatment_events")
    }
    expect_error({
        y <- add_continuous_data(x, timevar_data = list(A = points),
                                 treatment_events = c("A", "unknown"),
                                 covariate_events = NULL)
        prepare_rtmle_data(y)
    }, "treatment_events|unknown")
    expect_error({
        y <- add_continuous_data(x, timevar_data = list(A = points),
                                 covariate_events = NULL)
        prepare_rtmle_data(y)
    }, "treatment_events")
})

test_that("treatment values require explicit point measurements", {
    x <- continuous_treatment_input()
    for (bad in list(data.frame(id = 1:2, start_date = c(1, 0), end_date = c(2, 3)),
                     data.frame(id = 1:2, date = c(1, 2)))) {
        expect_error({
            y <- add_continuous_data(x, timevar_data = list(A = bad),
                                     treatment_events = "A", covariate_events = NULL)
            prepare_rtmle_data(y)
        }, "treatment|interval|value")
    }
})

test_that("multiple named opportunity streams are united and unchanged assignments retained", {
    x <- continuous_treatment_input()
    x <- add_continuous_data(
        x, outcome_data = data.frame(id = 1, date = 3),
        censored_data = data.frame(id = 2, date = 3, event = "tauend"),
        timevar_data = list(A = continuous_treatment_points(),
                            visits = data.frame(id = 1:2, date = 2)),
        treatment_events = c("A", "visits"), covariate_events = NULL)
    expect_identical(x$continuous_event_roles$treatment, c("A", "visits"))
    # Omitting the role on a subsequent call preserves it; event-only
    # covariates continue to be valid.
    x <- add_continuous_data(x, bleeding = data.frame(id = 1, date = 1),
                             covariate_events = "bleeding")
    expect_identical(x$continuous_event_roles$treatment, c("A", "visits"))
    x <- prepare_rtmle_data(x, min_events = 0)
    expect_equal(x$prepared_data$time_1, c(1, 1))
    expect_equal(x$prepared_data$time_2, c(2, 2))
    expect_equal(as.character(x$prepared_data$event_1), c("A", "A"))
    expect_equal(as.character(x$prepared_data$event_2), c("A", "A"))
    expect_equal(x$prepared_data$A_1, c(0, 1))
    expect_equal(x$prepared_data$A_2, c(1, 0))
    expect_equal(x$prepared_data$bleeding_1, c(1, 0))
    expect_equal(x$prepared_data$time_3, c(3, 3))
})

test_that("each declared opportunity needs an observed treatment value", {
    make_analysis <- function(A) {
        x <- continuous_treatment_input()
        add_continuous_data(
            x, outcome_data = data.frame(id = 1, date = 3),
            timevar_data = list(A = A,
                                visits = data.frame(id = 1, date = 1.5)),
            treatment_events = c("A", "visits"), covariate_events = NULL)
    }
    points <- continuous_treatment_points()
    expect_error(prepare_rtmle_data(make_analysis(points), min_events = 0),
                 "treatment.*(observ|value)|observ.*treatment")

    # An explicit unchanged value records the treatment decision at the visit.
    points <- rbind(points, data.frame(id = 1, date = 1.5, value = 0))
    x <- prepare_rtmle_data(make_analysis(points), min_events = 0)
    expect_equal(x$prepared_data[id == 1, time_2], 1.5)
    expect_equal(as.character(x$prepared_data[id == 1, event_2]), "A")
    expect_equal(x$prepared_data[id == 1, A_2], 0)
})
