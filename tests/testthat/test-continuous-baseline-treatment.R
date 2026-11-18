baseline_treatment_example <- function(continuous, baseline, A, L = data.frame(id = 1, date = 0, value = 10)) {
    x <- rtmle_init(c(0, 1, 3), name_id = "id", name_outcome = "Y",
                    continuous_time = continuous)
    x <- add_baseline_data(x, baseline)
    if (continuous) {
        x <- add_continuous_data(x, outcome_data = data.frame(id = 1, date = 3),
                                 timevar_data = list(A = A, L = L),
                                 treatment_events = "A", covariate_events = "L")
        prepare_rtmle_data(x, min_events = 0)
    } else {
        x <- add_long_data(x, outcome_data = data.frame(id = 1, date = 3),
                           timevar_data = list(A = A, L = L))
        x <- discretize_data(x, start_followup_date = 0,
                             mappings = list(A = "locf", L = "locf"), verbose = FALSE)
        prepare_rtmle_data(x)
    }
}

test_that("treatment at baseline carries forward the latest prebaseline value", {
    baseline <- data.frame(id = 1, age = 50)
    A <- data.frame(id = c(1, 1, 1), date = c(-3, -5, 2), value = c(1, 0, 0))
    discrete <- baseline_treatment_example(FALSE, baseline, A)
    continuous <- baseline_treatment_example(TRUE, baseline, A)
    expect_equal(discrete$prepared_data$A_0, 1)
    expect_equal(continuous$prepared_data$A_0, discrete$prepared_data$A_0)
    # Later treatment decisions must still be retained.
    expect_equal(continuous$prepared_data$time_1, 2)
    expect_equal(continuous$prepared_data$A_1, 0)
})

test_that("baseline treatment history is resolved separately for each subject", {
    make <- function(continuous) {
        x <- rtmle_init(c(0, 1, 3), name_id = "id", name_outcome = "Y",
                        name_censoring = "C", continuous_time = continuous)
        x <- add_baseline_data(x, data.frame(id = 1:2, age = c(50, 60)))
        data <- list(
            outcome_data = data.frame(id = 1, date = 3),
            censored_data = data.frame(id = 2, date = 3, event = "tauend"),
            timevar_data = list(
                A = data.frame(id = c(1, 1, 2), date = c(-3, 2, 2), value = c(1, 0, 1)),
                L = data.frame(id = 1:2, date = 0, value = c(10, 20))))
        if (continuous) {
            x <- do.call(add_continuous_data, c(list(x = x), data,
                                               list(treatment_events = "A", covariate_events = "L")))
            prepare_rtmle_data(x, min_events = 0)
        } else {
            x <- do.call(add_long_data, c(list(x = x), data))
            x <- discretize_data(x, start_followup_date = 0,
                                 mappings = list(A = "locf", L = "locf"), verbose = FALSE)
            prepare_rtmle_data(x)
        }
    }
    discrete <- make(FALSE)
    continuous <- make(TRUE)
    expect_equal(discrete$prepared_data$A_0, c(1, 0))
    expect_equal(continuous$prepared_data$A_0, discrete$prepared_data$A_0)
    expect_equal(continuous$prepared_data$A_1, c(0, 1))
})

test_that("treatment with no baseline history initializes to zero", {
    baseline <- data.frame(id = 1, age = 50)
    A <- data.frame(id = 1, date = 2, value = 1)
    discrete <- baseline_treatment_example(FALSE, baseline, A)
    continuous <- baseline_treatment_example(TRUE, baseline, A)
    expect_equal(discrete$prepared_data$A_0, 0)
    expect_equal(continuous$prepared_data$A_0, discrete$prepared_data$A_0)
    expect_equal(continuous$prepared_data$A_1, 1)
})

test_that("supplied baseline treatment takes precedence over a conflicting time-zero record", {
    A <- data.frame(id = c(1, 1), date = c(0, 2), value = c(0, 0))
    baseline <- data.frame(id = 1, age = 50, A_0 = 1)
    discrete <- baseline_treatment_example(FALSE, baseline, A)
    continuous <- baseline_treatment_example(TRUE, baseline, A)
    expect_equal(discrete$prepared_data$A_0, 1)
    expect_equal(continuous$prepared_data$A_0, discrete$prepared_data$A_0)

    # The user also accepts the unsuffixed name A as an initial-value alias.
    alias <- baseline_treatment_example(TRUE, data.frame(id = 1, age = 50, A = 1), A)
    expect_equal(alias$prepared_data$A_0, 1)
    expect_equal(alias$prepared_data$A_1, 0)
})

test_that("unsuffixed covariate values initialize L and conflicting baseline aliases are rejected", {
    baseline <- data.frame(id = 1, age = 50, L = 10)
    A <- data.frame(id = c(1, 1), date = c(0, 2), value = c(0, 1))
    # No L observation at time zero: it must come from baseline L.
    L <- data.frame(id = 1, date = 2, value = 12)
    continuous <- baseline_treatment_example(TRUE, baseline, A, L)
    expect_equal(continuous$prepared_data$L_0, 10)
    expect_equal(continuous$prepared_data$L_1, 12)
    baseline$L_0 <- 11
    expect_error(baseline_treatment_example(TRUE, baseline, A, L), "conflict")
})
