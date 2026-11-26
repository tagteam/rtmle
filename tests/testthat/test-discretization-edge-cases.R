test_that("discretization reports incomplete event histories and invalid mappings", {
    x <- rtmle_init(time_grid = 0:2, name_id = "id", name_outcome = "Y")
    expect_error(discretize_data(x), "need the baseline data")
    x <- add_baseline_data(x, data.frame(id = 1:2))
    expect_error(discretize_data(x), "long_data has length 0", fixed = TRUE)
    empty <- x
    empty$long_data <- list(timevar_data = list())
    expect_error(discretize_data(empty), "does not contain outcome data")
    x <- add_long_data(x, Z = data.frame(id = 1, date = 1, value = 3))
    expect_error(discretize_data(x), "does not contain outcome data")
    x <- add_long_data(x, outcome_data = data.frame(id = 1:2, date = c(1, 2)))

    expect_warning(discretize_data(x, breaks = 0:2), "obsolete")
    expect_error(discretize_data(x, 0, NULL, FALSE, 0, "locf"), "arguments must be named")
    expect_error(discretize_data(x, Z = 12), "character string or a named list")
    expect_error(discretize_data(x, W = "locf"), "longformat variable 'W'")
    expect_error(discretize_data(x, Z = "undefined"), "Unknown method for 'Z'")
    expect_error(discretize_data(x, Z = list(method = identity)), "required input columns")
    expect_error(discretize_data(x, Z = "exposure_time"), "start_date, end_date")
    expect_error(discretize_data(x, start_followup_date = "typo"), "name .* of a variable")
    expect_error(discretize_data(x, Z = list(method = "locf", baseline_lookback = -1)),
                 "Mapping 'Z' baseline_lookback")

    # Support custom functions supplied by either 'method' or 'fun'.
    double_measurement <- function(data, grid, name, id) {
        data <- data.table::copy(data)
        data$value <- 2 * data$value
        map_data_to_grid("measurement", data, grid, name, id)
    }
    mapped <- discretize_data(x, baseline_lookback = NULL,
                              by_method = list(variable = "Z", method = double_measurement,
                                               columns = c("date", "value")),
                              by_fun = list(variable = "Z", fun = double_measurement,
                                            columns = c("date", "value")))
    expect_equal(mapped$data$timevar_data$by_method$by_method_1, c(6, NA))
    expect_equal(mapped$data$timevar_data$by_fun$by_fun_1, c(6, NA))
})

test_that("invalid layouts and inconsistent date classes are identified", {
    x <- rtmle_init(time_grid = 0:2, name_id = "id", name_outcome = "Y")
    x <- add_baseline_data(x, data.frame(id = 1:2, start = as.Date("2020-01-01")))
    x <- add_long_data(x, outcome_data = data.frame(id = 1:2, date = c(1, 2)))
    expect_error(discretize_data(x, start_followup_date = "start"), "requires Date-valued")

    # Older or manually constructed objects may bypass add_long_data's checks.
    broken <- x
    broken$long_data$outcome_data <- data.frame(id = 1:2)
    expect_error(discretize_data(broken), "must contain a 'date' column")
    broken$long_data$timevar_data <- list(L = data.frame(id = 1, value = 5))
    expect_error(discretize_data(broken), "Could not determine method")
    expect_error(discretize_data(broken, L = list(method = identity, columns = "value")),
                 "must contain date or start_date/end_date")
    broken$long_data$timevar_data <- c(broken$long_data$timevar_data, broken$long_data$timevar_data)
    expect_error(discretize_data(broken), "Duplicated names")

    x_char <- add_long_data(x, outcome_data = data.frame(id = 1L, date = "2020-01-02"),
                            Z = data.frame(id = 1L, date = "2020-01-01"))
    expect_error(discretize_data(x_char), "must be numeric or Date")
    bad_interval <- data.frame(id = 1L, start_date = 0, end_date = "1")
    expect_error(discretize_data(add_long_data(x, A = bad_interval)), "must have the same class")
    bad_interval$start_date <- "0"
    expect_output(expect_error(discretize_data(add_long_data(x, A = bad_interval)),
                               "need to have the same class"), "class_of_date_variable")

    starts <- as.Date(c("2020-01-01", "2020-01-04"))
    x <- add_long_data(x, outcome_data = data.frame(id = 1:2, date = starts + 2))
    expect_error(discretize_data(add_baseline_data(x, data.frame(id = 1:2, start = 0)),
                                 start_followup_date = "start"), "Date-valued")
    expect_error(discretize_data(add_baseline_data(x, data.frame(id = 1:2, start = c(starts[1], NA))),
                                 start_followup_date = "start"), "cannot contain missing")
    x <- add_long_data(x, A = data.frame(id = 3L, date = starts[1]))
    expect_error(discretize_data(x, start_followup_date = "start"), "not present in the baseline")
    x$progress$substracted_start_followup_date <- TRUE
    expect_error(discretize_data(x, start_followup_date = "start"), "marked as converted")
})

test_that("insufficient follow-up and censoring before every outcome are flagged", {
    x <- rtmle_init(time_grid = 0:3, name_id = "id", name_outcome = "Y",
                    name_competing = "D", name_censoring = "C")
    x <- add_baseline_data(x, data.frame(id = 1:2))
    x <- add_long_data(x, outcome_data = data.frame(id = 1:2, date = 1))
    expect_error(discretize_data(x), "maximal followup time")
    x <- add_long_data(x, outcome_data = data.frame(id = 1:2, date = 6),
                       censored_data = data.frame(id = 1:2, date = 3))
    expect_error(discretize_data(x), "all outcome dates are administratively censored")
    # Person 2 has no recorded time of any event.
    incomplete <- add_long_data(x, outcome_data = data.frame(id = 1, date = 1),
                                censored_data = data.frame(id = 1, date = 3),
                                competing_data = data.frame(id = 1, date = 2))
    expect_error(discretize_data(incomplete), "Missing values in end of followup")
})

test_that("widen_outcome can read original event dates and rejects repeated terminal dates", {
    x <- rtmle_init(time_grid = 0:2, name_id = "id", name_outcome = "Y",
                    name_censoring = "C", name_competing = "D")
    expect_error(rtmle:::widen_outcome(x), "does not contain outcome data")
    x <- add_long_data(x, outcome_data = data.frame(id = 1, date = 1.5))
    expect_error(rtmle:::widen_outcome(x), "grid is required")
    grid <- data.table::data.table(id = rep(1:2, each = 3), interval = rep(0:2, 2),
                                   start_interval = c(0, 0, 1), end_interval = 0:2)
    duplicated_event <- data.table::data.table(id = c(2, 2), date = c(0.5, 1.5))
    expect_error(rtmle:::widen_outcome(x, grid = grid, censored_data = duplicated_event),
                 "Duplicated person id in dates for censoring")
    expect_error(rtmle:::widen_outcome(x, grid = grid, competing_data = duplicated_event),
                 "Duplicated person id in dates for competing")
    w <- rtmle:::widen_outcome(x, grid = grid, competing_data = duplicated_event[1L])
    expect_equal(w$D_1, c(0, 1))
    expect_equal(w$Y_2, c(1, 0))
})

test_that("chronic disease progresses from never to acute to chronic after an event", {
    x <- rtmle_init(time_grid = c(0, 10, 200), name_id = "id", name_outcome = "Y")
    onset <- as.Date("2020-01-01")
    x <- add_baseline_data(x, data.frame(id = 1:2, start = onset))
    x <- add_long_data(x, outcome_data = data.frame(id = 1:2, date = onset + 210),
                       bleeding = data.frame(id = 1, date = onset + 1),
                       other = data.frame(id = integer(), date = as.Date(character())))
    x <- discretize_data(x, start_followup_date = "start",
                         status = list(method = "chronic_disease", variable = "bleeding"),
                         other = list(method = "event", fun = NULL))
    w <- x$data$timevar_data$status
    expect_equal(as.character(w$status_0), c("never", "never"))
    expect_equal(as.character(w$status_1), c("acute", "never"))
    expect_equal(as.character(w$status_2), c("chronic", "never"))
    prepared <- prepare_rtmle_data(x)
    expect_false("start" %in% names(prepared$prepared_data))
    expect_false("start" %in% prepared$names$name_baseline_covariates)
})
