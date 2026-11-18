continuous_long_fixture <- function() {
    list(names = list(id = "person", time = "visit", outcome = "stroke",
                      competing = "death", censoring = "dropout"),
         time_grid_scale = c(0, 4),
         continuous_event_roles = list(treatment = "A", covariate = "L"),
         data = list(baseline_data = data.frame(person = 1:2, A = c(0, 1),
                                               L = c(10, 20), age = c(50, 60))),
         long_data = list(
             outcome_data = data.frame(person = 1, visit = 3),
             timevar_data = list(
                 A = data.frame(person = c(1, 1, 2), visit = c(1, 2, 1),
                                value = c(0, 1, 1)),
                 L = data.frame(person = c(1, 2), visit = c(1, 2), value = c(11, 21)))))
}

test_that("separate continuous streams retain opportunities and carry histories", {
    x <- continuous_long_fixture()
    before <- data.table::copy(x)
    out <- continuous_data_from_long(x)
    events <- out$data$timevarying_data
    expect_identical(x, before)
    expect_equal(out$time_covariates, c("A", "L"))
    expect_setequal(out$baseline_covariates, c("A_0", "L_0", "age"))
    expect_equal(out$data$baseline_data$A_0, c(0, 1))
    expect_equal(events[id == 1, time], c(1, 2, 3))
    expect_equal(events[id == 1, event], c("A", "A", "Y"))
    expect_equal(events[id == 1, A], c(0, 1, 1))
    expect_equal(events[id == 1, L], c(11, 11, 11))
    expect_equal(events[id == 2, time], c(1, 2, 5))
    expect_equal(events[id == 2, event], c("A", "L", "tauend"))
    expect_equal(events[id == 2, L], c(20, 21, 21))
    x$data$baseline_data$L <- c(FALSE, TRUE)
    x$long_data$timevar_data$L$value <- c(1, 0)
    expect_equal(continuous_data_from_long(x)$data$timevarying_data[id == 1, L],
                 c(1, 1, 1))
})

test_that("continuous terminals use the earliest event when times differ", {
    x <- continuous_long_fixture()
    x$long_data$outcome_data <- data.frame(person = 1:2, visit = c(3, 5))
    x$long_data$competing_data <- data.frame(person = 1:2, visit = c(4, 3.5))
    x$long_data$censored_data <- data.frame(person = 1:2, visit = c(5, 4.5))
    events <- continuous_data_from_long(x)$data$timevarying_data
    expect_equal(events[id == 1, event], c("A", "A", "Y"))
    expect_equal(events[id == 1, L], c(11, 11, 11))
    expect_equal(events[id == 2, event], c("A", "L", "D"))
    expect_equal(events[id == 2, time], c(1, 2, 3.5))
})

test_that("terminal events cannot tie other terminal, treatment or covariate events", {
    # Ties between outcome, competing and censoring records are rejected even
    # if an earlier terminal record would otherwise remove them from follow-up.
    x <- continuous_long_fixture()
    x$long_data$competing_data <- data.frame(person = 1, visit = 3)
    expect_error(continuous_data_from_long(x), "tied.*terminal")
    x$long_data$competing_data$visit <- 4
    x$long_data$censored_data <- data.frame(person = 1, visit = 4)
    expect_error(continuous_data_from_long(x), "tied.*terminal")

    # Treatment and separately declared covariate events cannot coincide with
    # any terminal event, including inferred administrative ends.
    x <- continuous_long_fixture()
    x$long_data$outcome_data$visit <- 2
    expect_error(continuous_data_from_long(x), "tied.*terminal")
    x <- continuous_long_fixture()
    x$long_data$competing_data <- data.frame(person = 2, visit = 2)
    expect_error(continuous_data_from_long(x), "tied.*terminal")
    x <- continuous_long_fixture()
    x$long_data$timevar_data$A <- rbind(x$long_data$timevar_data$A,
                                      data.frame(person = 2, visit = 5, value = 1))
    expect_error(continuous_data_from_long(x), "tied.*terminal")
})

test_that("terminal-time covariate measurements are not separate events", {
    x <- continuous_long_fixture()
    x$long_data$timevar_data$L <- rbind(x$long_data$timevar_data$L,
                                      data.frame(person = 1, visit = 3, value = 13))
    x$continuous_event_roles$covariate <- data.frame(person = 2, visit = 2)
    events <- continuous_data_from_long(x)$data$timevarying_data
    expect_equal(events[id == 1, event], c("A", "A", "Y"))
    expect_equal(events[id == 1, L], c(11, 11, 13))

    # Equal times across different people are not ties within a subject.
    x$long_data$competing_data <- data.frame(person = 2, visit = 3)
    expect_equal(continuous_data_from_long(x)$data$timevarying_data[id == 2, event],
                 c("A", "L", "D"))
})

test_that("one terminal event table row may carry covariate measurements", {
    x <- rtmle_init(c(0, 4), name_id = "id", name_time = "time",
                    name_outcome = "Y", name_censoring = "C", continuous_time = TRUE)
    x <- add_baseline_data(x, data.frame(id = 1:2, A = c(0, 1), L = c(10, 20)))
    events <- data.frame(id = c(1, 1, 2, 2), time = c(1, 3, 1, 3),
                         event = c("A", "Y", "A", "C"),
                         A = c(1, 1, 0, 0), L = c(11, 13, 21, 23))
    x <- add_continuous_data(x, event_data = events, treatment_events = "A",
                             covariate_events = NULL, time_covariates = "L")
    prepared <- prepare_rtmle_data(x, min_events = 0)$prepared_data
    expect_equal(as.character(prepared$event_2), c("Y", "C"))
    expect_equal(prepared$L_2, c(13, 23))
})

test_that("continuous administrative endpoints and later follow-up are preserved", {
    x <- continuous_long_fixture()
    x$long_data$outcome_data <- NULL
    x$long_data$censored_data <- data.frame(person = 1:2, visit = c(5, 6),
                                           event = c("C", "tauend"))
    events <- continuous_data_from_long(x)$data$timevarying_data
    expect_equal(events[id == 1, tail(event, 1)], "C")
    expect_equal(events[id == 1, tail(time, 1)], 5)
    expect_equal(events[id == 2, tail(event, 1)], "tauend")
    expect_equal(events[id == 2, tail(time, 1)], 6)
    x$long_data$censored_data$event <- NULL
    expect_true(all(continuous_data_from_long(x)$data$timevarying_data[
        time > 4, event] == "C"))
})

test_that("time-zero values initialize baseline without overwriting supplied values", {
    x <- continuous_long_fixture()
    x$data$baseline_data$A <- NULL
    x$long_data$timevar_data$A <- rbind(data.frame(person = 1:2, visit = 0, value = c(0, 1)),
                                      x$long_data$timevar_data$A)
    expect_equal(continuous_data_from_long(x)$data$baseline_data$A_0, c(0, 1))
    x$data$baseline_data$A_0 <- c(1, 1)
    expect_error(continuous_data_from_long(x), "conflict")
})

test_that("baseline-only treatment and event-only histories use existing data conventions", {
    x <- continuous_long_fixture()
    x$long_data$timevar_data$A <- data.frame(person = integer(), visit = numeric(), value = numeric())
    x$continuous_event_roles$treatment <- "A"
    x$long_data$timevar_data$L$value <- NULL
    x$data$baseline_data$L <- NULL
    out <- continuous_data_from_long(x)
    expect_setequal(out$time_covariates, c("L", "A"))
    expect_equal(out$data$baseline_data$L_0, c(0, 0))
    expect_equal(out$data$timevarying_data[id == 1, event], c("L", "Y"))
    expect_equal(out$data$timevarying_data[id == 1, A], c(0, 0))
    expect_equal(out$data$timevarying_data[id == 1, L], c(1, 1))
})

test_that("continuous stream errors are rejected before preparation", {
    x <- continuous_long_fixture()
    x$long_data$timevar_data$L$person[1] <- 99
    expect_error(continuous_data_from_long(x), "Unknown subject")
    x <- continuous_long_fixture()
    x$long_data$timevar_data$A <- rbind(x$long_data$timevar_data$A,
                                      x$long_data$timevar_data$A[1, ])
    expect_error(continuous_data_from_long(x), "Duplicated")
    x <- continuous_long_fixture()
    x$data$baseline_data$L <- NULL
    expect_error(continuous_data_from_long(x), "initial")
    x <- continuous_long_fixture()
    x$long_data$timevar_data$A$visit <- as.Date("2020-01-01") + 1:3
    expect_error(continuous_data_from_long(x), "numeric durations")
    x <- continuous_long_fixture()
    x$long_data$timevar_data$A <- data.frame(person = 1, start_date = 1, end_date = 2)
    expect_error(continuous_data_from_long(x), "treatment|interval")
    x <- continuous_long_fixture()
    x$long_data$censored_data <- data.frame(person = 1, visit = 3, event = "Y")
    expect_error(continuous_data_from_long(x), "C or tauend")
})

test_that("treatment measurements do not define treatment opportunities", {
    x <- continuous_long_fixture()
    x$long_data$timevar_data$visits <- data.frame(person = 1, visit = 2)
    x$continuous_event_roles <- list(
        treatment = "visits",
        covariate = data.frame(person = c(1, 2, 2), visit = c(1, 1, 2)))
    events <- continuous_data_from_long(x)$data$timevarying_data
    expect_equal(events[id == 1, event], c("L", "A", "Y"))
    expect_equal(events[id == 2, event], c("L", "L", "tauend"))
    x$long_data$timevar_data$A$value[1] <- 1
    expect_error(continuous_data_from_long(x), "declared treatment opportunity")
    x <- continuous_long_fixture()
    x$continuous_event_roles$covariate <- data.frame(person = 1, visit = 1)
    expect_error(continuous_data_from_long(x), "Undeclared measurement")
})

test_that("interval boundaries supply state without creating opportunities", {
    x <- continuous_long_fixture()
    x$long_data$timevar_data$L <- data.frame(person = 1, visit = 1.5, value = 99)
    attr(x$long_data$timevar_data$L, "continuous_interval") <- TRUE
    x$continuous_event_roles["covariate"] <- list(NULL)
    events <- continuous_data_from_long(x)$data$timevarying_data
    expect_equal(events[id == 1, time], c(1, 2, 3))
    expect_equal(events[id == 1, L], c(10, 99, 99))
    x$continuous_event_roles$covariate <- "L"
    expect_error(continuous_data_from_long(x), "explicit event times")
    x$long_data$timevar_data$A <- data.frame(person = 1, visit = 1.5, value = 1)
    attr(x$long_data$timevar_data$A, "continuous_interval") <- TRUE
    x$continuous_event_roles <- list(treatment = "A", covariate = NULL)
    expect_error(continuous_data_from_long(x), "treatment|interval|point")
})

test_that("continuous event roles must be declared explicitly", {
    x <- continuous_long_fixture()
    x$continuous_event_roles <- NULL
    expect_error(continuous_data_from_long(x), "Declare treatment_events and covariate_events")
})
