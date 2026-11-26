test_that("a linear cumulative hazard gives the inverse exponential transform", {
    # H(t) = 2t, with administrative follow-up ending at time 4.
    n <- 16
    ratios <- rep(c(0.1, 1, 2, 4), length.out = n)
    set.seed(408)
    exponential_draw <- stats::rexp(n)
    set.seed(408)
    observed <- rtmle:::reventtime(n, breaks = c(1, 2, 4),
                                   cumhazard = c(2, 4, 8), hazardratio = ratios,
                                   entrytime = 1.5, decimals = 3)
    expect_equal(observed, round(pmin(4, 1.5 + exponential_draw / (2 * ratios)), 3))
    expect_error(rtmle:::reventtime(breaks = 0:2, cumhazard = 0:2, hazardratio = 1),
                 "specify sample size")
    expect_error(rtmle:::reventtime(2, 0:2, 0:2, hazardratio = c(1, 0)),
                 "cannot contain zeros")
})

test_that("cohort simulation respects scheduled visits, hooks and fixed interventions", {
    cohort <- simulate_cohort(
        n = 2, seed = 48, max_follow = 4,
        baseline_variables = data.table::data.table(age = c(50, 60)),
        post_baseline_variables_hook = function(X) {
            X$visit_count <- 0L
            X
        },
        baseline_visit = list(A = "constant"),
        post_baseline_visit_hook = function(X) {
            X$A <- 0L
            X
        },
        visit_schedule = list(mean = 8, sd = 0, skip = 0, schedule = c(1, 5)),
        next_visit_hook = function(next_visit, last_entry) next_visit / 4,
        visit_events = list(A = "binomial"),
        intermediate_events = list(bleeding = "constant"),
        pre_intermediate_events_hook = function(intermediate_events_model, event_history) {
            lava::distribution(intermediate_events_model, "bleeding") <- lava::constant.lvm(100)
            intermediate_events_model
        },
        post_intermediate_events_hook = function(update_information, event_history) {
            update_information$visit_count <- update_information$visit_count + 1L
            update_information
        },
        post_visit_hook = function(update_information, event_history) {
            update_information$visit_score <- update_information$visit_count * 10L
            update_information
        },
        absorbing_events = list(death = "constant", dropout = "constant"),
        absorbing_events_hook = function(absorbing_event_model, event_history) {
            if (all(event_history$visit_count == 2)) {
                lava::distribution(absorbing_event_model, "death") <- lava::constant.lvm(0.5)
            }
            absorbing_event_model
        },
        parameter_values = list(intercept_A = 0, intercept_bleeding = 0.25,
                                intercept_death = 50, intercept_dropout = 0.1),
        intervention = c(A = 1)
    )
    expect_equal(cohort$time, rep(c(0, 1, 3, 3.5), 2))
    expect_equal(as.character(cohort$event), rep(c("baseline", "visit", "visit", "death"), 2))
    visits <- cohort[as.character(cohort$event) == "visit", ]
    expect_equal(visits$A, rep(1, 4))
    expect_equal(visits$visit_score, rep(c(10L, 20L), 2))
    # Intervention names must refer to visit treatments.
    expect_error(simulate_cohort(n = 2, visit_events = list(A = "binomial"),
                                 intervention = c(B = 1)), "names must match")
    expect_error(simulate_cohort(n = 2, visit_events = list(A = "binomial"),
                                 intervention = 1), "must be named")
})

test_that("register conversion merges consecutive exposed visits", {
    cohort <- simulate_cohort(n = 2, seed = 8, max_follow = 3,
                              visit_schedule = list(mean = 1, sd = 0, skip = 0),
                              visit_events = list(A = "constant"),
                              absorbing_events = list(death = "constant"),
                              parameter_values = list(intercept_A = 1, intercept_death = 10))
    reg <- register_format(cohort, treatment_variables = "A")
    expect_equal(reg$timevar_data$A$id, 1:2)
    expect_equal(reg$timevar_data$A$start_date, c(1, 1))
    expect_equal(reg$timevar_data$A$end_date, c(3, 3))
    expect_error(register_format(cohort, treatment_variables = NA_character_), "character vector")
    expect_error(register_format(cohort, treatment_variables = "unknown"), "not visit variables")
    expect_error(register_format(data.frame(id = 1)), "simulated_cohort")
})

test_that("long-data simulation assigns every visit treatment under intervention", {
    set.seed(182)
    population <- simulate_long_data(
        n = 30, number_visits = 4, interventional_distribution = TRUE,
        baseline_rate = list(L = 1:3650 / 1e6, D = 1:3650 / 1e6, Y = 1:3650 / 1000)
    )
    expect_equal(unique(population$A_0), 1)
    visits <- population[as.character(population$event) == "V", ]
    expect_gt(nrow(visits), 0L)
    expect_equal(unique(visits$A), 1)
    observed_censoring <- population[as.character(population$terminal_event) == "C", ]
    expect_gt(nrow(observed_censoring), 0L)
    expect_true(all(observed_censoring$terminal_time == 3650))
})

test_that("patients with no scheduled visits progress directly to an absorbing event", {
    cohort <- simulate_cohort(
        n = 3, seed = 42, max_follow = 20,
        absorbing_events = list(death = "constant"), parameter_values = list(intercept_death = 7)
    )
    expect_equal(cohort$time, rep(c(0, 7), 3))
    expect_equal(as.character(cohort$event), rep(c("baseline", "death"), 3))
    expect_equal(cohort$id, rep(1:3, each = 2))
})
