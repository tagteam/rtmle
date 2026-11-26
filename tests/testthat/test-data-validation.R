test_that("invalid baseline, long and wide inputs identify the offending field", {
    x <- rtmle_init(time_grid = 0:2, name_id = "pid", name_outcome = "Y")
    expect_error(add_baseline_data(x, 1), "data.frame")
    expect_error(add_baseline_data(x, data.frame(id = 1)), "pid")
    expect_error(add_baseline_data(x, data.frame(pid = c(1, 1))), "Duplicated")

    expect_error(add_long_data(x, timevar_data = data.frame(pid = 1, date = 1)),
                 "named list")
    expect_error(add_long_data(x, A = 1), "must be a data.frame")
    expect_error(add_long_data(x, outcome_data = 1), "outcome_data must be a data.frame")
    expect_error(add_long_data(x, NULL, NULL, NULL, NULL, data.frame(pid = 1)),
                 "arguments must be named")
    a <- data.frame(pid = 1, date = 0, value = 10)
    x <- add_long_data(x, A = a)
    a$value <- 25
    x <- add_long_data(x, A = a)
    expect_equal(x$long_data$timevar_data$A$value, 25)

    expect_error(add_wide_data(x, outcome_data = 1), "outcome_data must be a data.frame")
    expect_error(add_wide_data(x, outcome_data = data.frame(Y_1 = 0)), "pid")
    expect_error(add_wide_data(x, outcome_data = data.frame(pid = 1)), "Cannot find the outcome")
    expect_error(add_wide_data(x, outcome_data = data.frame(pid = 1:2, Y_0 = 0:1, Y_1 = 0)),
                 "variation in the time zero variable Y_0")
    expect_error(add_wide_data(x, outcome_data = data.frame(pid = 1:3, Y_1 = 0:2)), "must be binary")
    expect_error(add_wide_data(x, timevar_data = 1), "named list")
    expect_error(add_wide_data(x, timevar_data = data.frame(A_0 = 0)), "pid")
    expect_error(add_wide_data(x, timevar_data = list(A = data.frame(A_0 = 0))), "pid")
})

test_that("data preparation checks for missing outcome, competing and censoring nodes", {
    x <- rtmle_init(time_grid = 0:3, name_id = "id", name_outcome = "Y",
                    name_competing = "D", name_censoring = "C")
    expect_error(prepare_rtmle_data(x), "contains no subjects")
    wide <- data.frame(id = 1:2, Y_1 = 0, Y_2 = 0, Y_3 = 0,
                       C_1 = "uncensored", C_2 = "uncensored", C_3 = "uncensored",
                       D_1 = 0, D_2 = 0)
    for (field in c("Y_2", "C_2", "D_2")) {
        partial <- add_wide_data(x, outcome_data = wide[setdiff(names(wide), field)])
        expect_error(prepare_rtmle_data(partial), field)
    }
    # Objects supplied directly in wide format also undergo a time-zero check.
    x$data$outcome_data <- transform(wide, Y_0 = 0:1)
    expect_error(prepare_rtmle_data(x), "Y_0 which is not constant")
})

test_that("an observed outcome takes priority over censoring in the same interval", {
    x <- rtmle_init(time_grid = 0:1, name_id = "id", name_outcome = "Y",
                    name_censoring = "C")
    x <- add_wide_data(x, outcome_data = data.frame(
        id = 1:3, C_1 = c("censored", "censored", "uncensored"), Y_1 = c(1, 0, 0)))
    fit_data <- prepare_rtmle_data(x)$prepared_data
    expect_equal(as.character(fit_data$C_1), c("uncensored", "censored", "uncensored"))
    expect_equal(fit_data$Y_1, c(1, NA, 0))
})

test_that("fast casting either aggregates duplicates or identifies them", {
    tbl <- data.table::data.table(pid = c(2, 2, 1), interval = c(1, 1, 0), value = c(4, 8, 2))
    expect_error(fast_cast(tbl, id = "pid", name = "SBP"), "Processing variable SBP")
    expect_error(fast_cast(tbl, id = "pid"), "Duplicate")
    expect_error(fast_cast(tbl, id = "patient"), "id")
    expect_error(fast_cast(tbl, id = "pid", value_col = "score"), "value_col")
    expect_error(fast_cast(tbl[, c("pid", "value")], id = "pid"), "interval")
    wide <- fast_cast(tbl, id = "pid", fun_aggregate = mean)
    expect_equal(wide$pid, c(2, 1))
    expect_equal(wide[["0"]], c(NA, 2))
    expect_equal(wide[["1"]], c(6, NA))
})

test_that("event histories and exposure thresholds map to hand-calculated intervals", {
    grid <- data.table::data.table(id = 1L, interval = 0:3,
                                   start_interval = c(0, 0, 2, 4), end_interval = c(0, 2, 4, 6))
    bleeding <- data.table::data.table(id = 1L, date = 1)
    exposure <- data.table::data.table(id = 1L, start_date = 0, end_date = 3)
    map <- function(method, data, ...) {
        out <- map_data_to_grid(method, data, grid, name = "A", id = "id", ...)
        unlist(out[, -1, with = FALSE], use.names = FALSE)
    }
    expect_null(map_data_to_grid("event", NULL, grid, name = "A", id = "id"))
    expect_equal(map("time_since_event", bleeding, values = NULL), c(NA, 1, 3, 5))
    expect_equal(map("event_interval", bleeding), c(0, 1, 0, 0))
    expect_equal(map("exposure_time", exposure), c(1, 2, 1, 0))
    expect_equal(map("any_exposure", exposure), c(1, 1, 1, 0))
    expect_equal(map("has_exposure", exposure, threshold = 0.5), c(1, 1, 0, 0))
    expect_error(map("has_exposure", exposure, threshold = 1.1), "in \\[0, 1")
    expect_error(map("any_exposure", exposure, threshold = NA_real_), "single numeric")
    expect_error(map("locf", bleeding), "value")
})

test_that("one wide table and separate named covariate tables give the same histories", {
    x <- rtmle_init(time_grid = 0:2, name_id = "pid", name_outcome = "Y")
    outcomes <- data.frame(pid = 1:3, Y_1 = c(0, 1, 0), Y_2 = c(1, 1, 0))
    a <- data.frame(pid = 1:3, A_0 = 0:2, A_1 = c(2, 1, 0))
    l <- data.frame(pid = 1:3, L_0 = 3:1, L_1 = 4:2)
    combined <- add_wide_data(x, outcome_data = outcomes,
                              timevar_data = merge(a, l, by = "pid"))
    separate <- add_wide_data(x, outcome_data = outcomes, timevar_data = list(A = a, L = l))
    ref <- prepare_rtmle_data(separate)
    expect_equal(prepare_rtmle_data(combined)$prepared_data, ref$prepared_data)
    expect_equal(ref$prepared_data$A_1, c(2, NA, 0))
    expect_equal(ref$followup$last_interval, c(1, 0, 1))
    expect_warning(prepare_rtmle_data(ref), "Overwriting existing prepared data")
})
