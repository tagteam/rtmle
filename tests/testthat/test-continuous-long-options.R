continuous_options_object <- function(baseline, streams, ...) {
    x <- rtmle_init(c(0, 10), name_id = "subject", name_time = "visit",
                    name_outcome = "Y", continuous_time = TRUE)
    x <- add_baseline_data(x, baseline)
    add_long_data(x, timevar_data = streams, ...)
}

test_that("continuous histories use numeric time since baseline", {
    x <- continuous_options_object(
        data.frame(subject = 1:2),
        list(A = data.frame(subject = 1:2, visit = c(0, 3), value = c(0, 1))),
        outcome_data = data.frame(subject = 1, visit = 5),
        censored_data = data.frame(subject = 2, visit = 10))
    before <- data.table::copy(x)
    result <- continuous_normalize_long_data(x)
    expect_equal(result$long_data$timevar_data$A[visit == 0, subject], 1:2)
    expect_equal(result$long_data$timevar_data$A[visit > 0, visit], 3)
    expect_equal(result$long_data$outcome_data$visit, 5)
    expect_equal(result$long_data$censored_data$visit, 10)
    expect_identical(x, before)

    x$long_data$outcome_data[, visit := as.Date("2026-01-05")]
    expect_error(continuous_normalize_long_data(x), "numeric.*since baseline")
})

test_that("interval histories use the union of overlapping and adjacent intervals", {
    x <- continuous_options_object(
        data.frame(subject = 1:3),
        list(B = data.frame(subject = c(1, 1, 1, 1, 1, 2),
                            start_date = c(-3, 1, 2, 6, 9, 4),
                            end_date = c(2, 4, 6, 8, Inf, 7))))
    result <- continuous_normalize_long_data(x)$long_data$timevar_data$B
    expect_equal(result[subject == 1, visit], c(0, 8, 9))
    expect_equal(result[subject == 1, value], c(1, 0, 1))
    expect_equal(result[subject == 2, visit], c(0, 4, 7))
    expect_equal(result[subject == 2, value], c(0, 1, 0))
    expect_equal(result[subject == 3, value], 0)
    expect_false(anyDuplicated(result[, .(subject, visit)]) > 0)
    bad <- data.table::copy(x)
    bad$long_data$timevar_data$B[1, end_date := start_date]
    expect_error(continuous_normalize_long_data(bad), "end|interval")
    bad$long_data$timevar_data$B[1, end_date := NA_real_]
    expect_error(continuous_normalize_long_data(bad), "missing|NA")
})

test_that("point measurements at baseline and later are retained", {
    x <- continuous_options_object(
        data.frame(subject = 1:2),
        list(B = data.frame(subject = c(1, 1, 1, 1, 2, 2),
                            visit = c(-20, -5, -1, 3, -7, 4),
                            value = c(0, 1, 0, 1, 1, 0))))
    result <- continuous_normalize_long_data(x)$long_data$timevar_data$B
    expect_true(all(result$visit > 0))
    x$long_data$timevar_data$B[3, visit := 0]
    result <- continuous_normalize_long_data(x)$long_data$timevar_data$B
    expect_equal(result[subject == 1, visit], c(0, 3))
    expect_equal(result[subject == 1, value], c(0, 1))
    duplicate <- data.table::copy(x)
    duplicate$long_data$timevar_data$B <- data.table::rbindlist(list(
        duplicate$long_data$timevar_data$B,
        data.table::data.table(subject = 1, visit = 3, value = 0)))
    expect_error(continuous_normalize_long_data(duplicate), "conflict|duplicate")
})

test_that("event-only histories initialize and increment event counts", {
    x <- continuous_options_object(data.frame(subject = 1:3),
                                    list(B = data.frame(subject = c(1, 1, 2), visit = c(-2, 4, 5))))
    result <- continuous_normalize_long_data(x)$long_data$timevar_data$B
    expect_equal(result[subject == 1, visit], c(0, 4))
    expect_equal(result[subject == 1, value], c(0, 1))
    expect_equal(result[subject == 2, value], c(0, 1))
    expect_equal(result[subject == 3, value], 0)
})

test_that("interval status is right continuous at zero and retains source metadata", {
    x <- continuous_options_object(data.frame(subject = 1:3),
                                    list(B = data.frame(subject = 1:3,
                                                        start_date = c(-2, 0, -2),
                                                        end_date = c(0, 2, Inf))))
    result <- continuous_normalize_long_data(x)
    points <- result$long_data$timevar_data$B
    expect_equal(points[visit == 0, value], c(0, 1, 1))
    expect_true(isTRUE(attr(points, "continuous_interval")))
    expect_identical(result$continuous_interval_variables, "B")
    x$continuous_event_roles <- list(treatment = "A", covariate = "B")
    expect_error(continuous_normalize_long_data(x), "point|interval")
})

test_that("separate event roles use numeric times since baseline", {
    x <- continuous_options_object(data.frame(subject = 1:2),
                                    list(A = data.frame(subject = 1:2, visit = 0,
                                                        value = c(0, 1))))
    x$continuous_event_roles <- list(
        treatment = "A",
        covariate = data.frame(subject = c(1, 1), date = c(-1, 3)))
    original <- data.table::copy(x)
    result <- continuous_normalize_long_data(x)
    expect_identical(result$continuous_event_roles$treatment, "A")
    expect_equal(result$continuous_event_roles$covariate$visit, 3)
    expect_identical(x, original)
})

test_that("event counts respect supplied initial values", {
    x <- continuous_options_object(data.frame(subject = 1:2, B_0 = factor(c(1, 0))),
                                    list(B = data.frame(subject = 2, visit = 5)))
    result <- continuous_normalize_long_data(x)$long_data$timevar_data$B
    expect_equal(result[visit == 0, value], c(1, 0))
    expect_equal(result[subject == 2 & visit == 5, value], 1)
})

test_that("event counts include time-zero events and deduplicate repeated records", {
    events <- data.frame(
        subject = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4),
        visit = c(-4, -1, -1, 0, 2, 2, 5, -5, 1, 3, 7, -3, 0, 1))
    x <- continuous_options_object(data.frame(subject = 1:4), list(B = events))
    original <- data.table::copy(x)
    points <- continuous_normalize_long_data(x)$long_data$timevar_data$B
    expect_equal(points[subject == 1, visit], c(0, 2, 5))
    expect_equal(points[subject == 1, value], c(1, 2, 3))
    expect_equal(points[subject == 2, visit], c(0, 1, 3, 7))
    expect_equal(points[subject == 2, value], c(0, 1, 2, 3))
    expect_equal(points[subject == 3, value], 0)
    expect_equal(points[subject == 4, visit], c(0, 1))
    expect_equal(points[subject == 4, value], c(1, 2))
    expect_identical(x, original)
})

test_that("supplied baseline counts are not double-counted at time zero", {
    events <- data.frame(subject = c(1, 1, 1, 1, 2, 2),
                         visit = c(-2, 0, 1, 2, -1, 1))
    x <- continuous_options_object(data.frame(subject = 1:3, B_0 = c(5, 2, 0)),
                                    list(B = events))
    points <- continuous_normalize_long_data(x)$long_data$timevar_data$B
    expect_equal(points[subject == 1, value], c(5, 6, 7))
    expect_equal(points[subject == 2, value], c(2, 3))
    expect_equal(points[subject == 3, value], 0)
    too_small <- data.table::copy(x)
    too_small$data$baseline_data[subject == 1, B_0 := 0]
    expect_error(continuous_normalize_long_data(too_small), "baseline|count")
})

test_that("initial event counts must be finite nonnegative integers", {
    events <- data.frame(subject = 1, visit = 1)
    for (invalid in list(-1, 1.5, Inf, NA_real_, "wrong", factor("wrong"))) {
        x <- continuous_options_object(data.frame(subject = 1, B_0 = invalid), list(B = events))
        expect_error(continuous_normalize_long_data(x), "count|initial|baseline")
    }
})

test_that("explicit point values remain measurements rather than counts", {
    x <- continuous_options_object(data.frame(subject = 1, B_0 = 10),
                                    list(B = data.frame(subject = c(1, 1, 1),
                                                        visit = c(0, 2, 3), value = c(10, 14, 12))))
    points <- continuous_normalize_long_data(x)$long_data$timevar_data$B
    expect_equal(points$value, c(10, 14, 12))
})

test_that("identical point duplicates collapse but missing times are rejected", {
    x <- continuous_options_object(data.frame(subject = 1),
                                    list(B = data.frame(subject = c(1, 1), visit = c(2, 2), value = c(1, 1))))
    result <- continuous_normalize_long_data(x)$long_data$timevar_data$B
    expect_equal(nrow(result), 1L)
    x$long_data$timevar_data$B[1, visit := NA_real_]
    expect_error(continuous_normalize_long_data(x), "Missing|NA")
})
