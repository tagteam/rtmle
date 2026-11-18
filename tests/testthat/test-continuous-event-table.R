# Compare the single-table interface with the established separate streams.
event_table_fixture <- function() {
    baseline <- data.frame(id = 1:4, A = c(0, 1, 0, 1),
                           L1 = c(10, 20, 30, 40), L2 = c(100, 200, 300, 400),
                           age = c(40, 50, 60, 70))
    events <- data.frame(
        id = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4),
        time = c(1, 2, 3, 4, 1, 2, 4, 1, 4, 1, 2, 5),
        event = c("A", "L", "review", "Y", "visit", "A", "D", "L", "C", "review", "visit", "tauend"),
        A = c(1, NA, 1, NA, NA, 0, NA, NA, NA, 1, NA, NA),
        L1 = c(NA, 11, NA, NA, 22, NA, NA, 33, NA, NA, 44, NA),
        L2 = c(NA, 101, NA, NA, 202, NA, NA, 303, NA, NA, 404, NA)
    )
    list(baseline = baseline, events = events)
}

event_table_init <- function(baseline) {
    x <- rtmle_init(c(0, 5), name_id = "id", name_time = "time",
                    name_outcome = "Y", name_competing = "D", name_censoring = "C",
                    continuous_time = TRUE)
    add_baseline_data(x, baseline)
}

# Use A for the combined opportunities and L1/L2 for the shared measurement
# times, avoiding any extra variables in the established representation.
event_table_separate <- function(baseline, events) {
    points <- function(rows, variable) {
        data.frame(id = events$id[rows], time = events$time[rows], value = events[[variable]][rows])
    }
    treatment <- events$event %in% c("A", "review")
    measurement <- events$event %in% c("L", "visit")
    x <- event_table_init(baseline)
    add_continuous_data(x,
        outcome_data = events[events$event == "Y", c("id", "time")],
        competing_data = events[events$event == "D", c("id", "time")],
        censored_data = events[events$event %in% c("C", "tauend"), c("id", "time", "event")],
        timevar_data = list(A = points(treatment, "A"), L1 = points(measurement, "L1"),
                            L2 = points(measurement, "L2")),
        treatment_events = "A", covariate_events = c("L1", "L2"))
}

event_table_combined <- function(baseline, events, ...) {
    add_continuous_data(event_table_init(baseline), event_data = events,
                        treatment_events = c("A", "review"),
                        covariate_events = c("L", "visit"), ...)
}

test_that("event tables preserve opportunities, multiple measurements, and terminal events", {
    fixture <- event_table_fixture()
    original <- fixture$events
    combined <- event_table_combined(fixture$baseline, fixture$events,
                                     time_covariates = c("L1", "L2"))
    expect_identical(fixture$events, original)
    single <- prepare_rtmle_data(combined, min_events = 0)
    separate <- prepare_rtmle_data(event_table_separate(fixture$baseline, fixture$events), min_events = 0)

    # Compare all prepared values as well as checking the expected histories.
    a <- as.data.frame(single$prepared_data)
    b <- as.data.frame(separate$prepared_data)
    expect_setequal(names(a), names(b))
    expect_equal(a[, sort(names(a))], b[, sort(names(b))])
    expect_equal(single$continuous_prepared$info, separate$continuous_prepared$info)
    expect_equal(a$A_0, c(0, 1, 0, 1))
    expect_equal(a$L1_0, c(10, 20, 30, 40))
    expect_equal(a$L2_0, c(100, 200, 300, 400))
    expect_equal(as.character(a$event_1), c("A", "L", "L", "A"))
    expect_equal(a$A_1, c(1, 1, 0, 1))
    expect_equal(a$L1_1, c(10, 22, 33, 40))
    expect_equal(a$L2_1, c(100, 202, 303, 400))
    expect_equal(as.character(a$event_2), c("L", "A", "C", "L"))
    expect_equal(a$A_2, c(1, 0, 0, 1))
    expect_equal(a$L1_2, c(11, 22, 33, 44))
    expect_equal(a$L2_2, c(101, 202, 303, 404))
    expect_equal(as.character(a$event_3), c("A", "D", NA, "tauend"))
    expect_equal(as.character(a$event_4[1]), "Y")
    # A repeated assignment at review is still a treatment opportunity.
    expect_equal(a$A_3[1], 1)
})

test_that("covariate columns are inferred and A on non-treatment rows is ignored", {
    fixture <- event_table_fixture()
    x <- event_table_combined(fixture$baseline, fixture$events)
    explicit <- event_table_combined(fixture$baseline, fixture$events,
                                     time_covariates = c("L1", "L2"))
    inferred <- prepare_rtmle_data(x, min_events = 0)$prepared_data
    expected <- prepare_rtmle_data(explicit, min_events = 0)$prepared_data
    expect_equal(inferred, expected)

    # Treatment values on terminal rows do not update treatment.
    events <- fixture$events
    terminal <- events$event %in% c("Y", "D", "C", "tauend")
    events$A[terminal] <- 123
    changed <- event_table_combined(fixture$baseline, events, time_covariates = c("L1", "L2"))
    expect_equal(prepare_rtmle_data(changed, min_events = 0)$prepared_data, expected)
})

test_that("covariates can change at treatment events without a separate L event", {
    fixture <- event_table_fixture()
    events <- fixture$events
    treatment <- events$event %in% c("A", "review")
    events$L1[treatment] <- c(12, 13, NA, 45)
    events$L2[treatment] <- c(102, NA, 203, NA)

    # Only L and visit are declared as separate covariate events. An observed
    # value on an A or review row still updates that covariate at the A event.
    combined <- event_table_combined(fixture$baseline, events,
                                     time_covariates = c("L1", "L2"))
    prepared <- prepare_rtmle_data(combined, min_events = 0)
    a <- prepared$prepared_data
    expect_equal(as.character(a$event_1), c("A", "L", "L", "A"))
    expect_equal(a$L1_1, c(12, 22, 33, 45))
    expect_equal(a$L2_1, c(102, 202, 303, 400))
    expect_equal(a$L1_2[2], 22) # No new L1 measurement at this treatment visit.
    expect_equal(a$L2_2[2], 203)
    expect_equal(a$L1_3[1], 13)
    expect_equal(a$L2_3[1], 101) # The preceding L event supplied this value.

    # Treatment models and dynamic rules use the current, pre-decision value.
    # Subject 1 crosses the rule's threshold at the treatment visit.
    first_treatment <- a[which(as.character(a$event_1) == "A")]
    history <- rtmle:::continuous_history(
        first_treatment, prepared$continuous_prepared$time_covariates,
        prepared$continuous_prepared$baseline_covariates, "propensity", k = 1)
    expect_true("L1_1" %in% history)
    expect_equal(first_treatment$L1_1, c(12, 45))
    rule <- function(data, intervention_table, time_node) {
        data[[paste0("A_", time_node)]] <- as.numeric(data[[paste0("L1_", time_node)]] >= 12)
        data
    }
    analysis <- rtmle_init(c(0, 5), name_id = "id", name_time = "time",
                           name_outcome = "Y", continuous_time = TRUE)
    analysis <- regime(analysis, "rule", intervention = 1, intervene_function = rule)
    dynamic <- analysis$regimes$rule
    dynamic$continuous_baseline_covariates <- prepared$continuous_prepared$baseline_covariates
    dynamic$continuous_time_covariates <- prepared$continuous_prepared$time_covariates
    decision <- rtmle:::continuous_intervention_policy(first_treatment, dynamic, 1L)
    expect_equal(first_treatment$L1_0[1], 10)
    expect_equal(decision$value, c(1, 1))
    expect_false(any(decision$natural))

    # The separate-table interface already accepts point updates for each
    # covariate at treatment times, including partial measurements.
    points <- function(rows, variable) {
        data.frame(id = events$id[rows], time = events$time[rows], value = events[[variable]][rows])
    }
    measured <- events$event %in% c("A", "review", "L", "visit")
    separate <- add_continuous_data(event_table_init(fixture$baseline),
        outcome_data = events[events$event == "Y", c("id", "time")],
        competing_data = events[events$event == "D", c("id", "time")],
        censored_data = events[events$event %in% c("C", "tauend"), c("id", "time", "event")],
        timevar_data = list(
            A = points(treatment, "A"),
            L1 = points(measured & !is.na(events$L1), "L1"),
            L2 = points(measured & !is.na(events$L2), "L2")),
        treatment_events = "A", covariate_events = c("L1", "L2"))
    b <- prepare_rtmle_data(separate, min_events = 0)$prepared_data
    expect_equal(a, b)
})

test_that("nonfinite covariate values at treatment events are rejected", {
    fixture <- event_table_fixture()
    for (value in c(Inf, -Inf, NaN)) {
        events <- fixture$events
        events$L1[1] <- value
        expect_error(event_table_combined(fixture$baseline, events,
                                          time_covariates = c("L1", "L2")),
                     "Nonfinite.*L1")
    }
})

test_that("covariates at terminal events are retained without extra L events", {
    fixture <- event_table_fixture()
    events <- fixture$events
    terminal <- events$event %in% c("Y", "D", "C", "tauend")
    events$L1[terminal] <- c(14, 24, 34, 45)
    events$L2[terminal] <- c(104, 204, 304, 405)
    combined <- event_table_combined(fixture$baseline, events,
                                     time_covariates = c("L1", "L2"))
    a <- prepare_rtmle_data(combined, min_events = 0)$prepared_data
    expect_equal(as.character(a$event_4[1]), "Y")
    expect_equal(c(a$L1_4[1], a$L2_4[1]), c(14, 104))
    expect_equal(as.character(a$event_3[2]), "D")
    expect_equal(c(a$L1_3[2], a$L2_3[2]), c(24, 204))
    expect_equal(as.character(a$event_2[3]), "C")
    expect_equal(c(a$L1_2[3], a$L2_2[3]), c(34, 304))
    expect_equal(as.character(a$event_3[4]), "tauend")
    expect_equal(c(a$L1_3[4], a$L2_3[4]), c(45, 405))

    # The same terminal-time values can be supplied as separate point streams,
    # without declaring those terminal times as covariate events.
    points <- function(rows, variable) {
        data.frame(id = events$id[rows], time = events$time[rows], value = events[[variable]][rows])
    }
    treatment <- events$event %in% c("A", "review")
    measurement <- events$event %in% c("L", "visit")
    separate <- add_continuous_data(event_table_init(fixture$baseline),
        outcome_data = events[events$event == "Y", c("id", "time")],
        competing_data = events[events$event == "D", c("id", "time")],
        censored_data = events[events$event %in% c("C", "tauend"), c("id", "time", "event")],
        timevar_data = list(
            A = points(treatment, "A"),
            L1 = points(measurement | terminal, "L1"),
            L2 = points(measurement | terminal, "L2")),
        treatment_events = "A", covariate_events = events[measurement, c("id", "time")])
    b <- prepare_rtmle_data(separate, min_events = 0)$prepared_data
    expect_equal(a, b)
})

test_that("nonfinite covariate values at terminal events are rejected", {
    fixture <- event_table_fixture()
    for (value in c(Inf, -Inf, NaN)) {
        events <- fixture$events
        events$L1[events$event == "Y"] <- value
        expect_error(event_table_combined(fixture$baseline, events,
                                          time_covariates = c("L1", "L2")),
                     "Nonfinite.*L1")
    }
})

test_that("one event label can record a treatment and covariate update together", {
    fixture <- event_table_fixture()
    review <- fixture$events$event == "review"
    fixture$events$L1[review] <- c(12, 41)
    fixture$events$L2[review] <- c(102, 401)
    x <- add_continuous_data(event_table_init(fixture$baseline), event_data = fixture$events,
                             treatment_events = c("A", "review"),
                             covariate_events = c("L", "visit", "review"),
                             time_covariates = c("L1", "L2"))
    wide <- prepare_rtmle_data(x, min_events = 0)$prepared_data
    # A simultaneous measurement is available at the treatment event.
    expect_equal(as.character(wide$event_3[1]), "A")
    expect_equal(wide$L1_3[1], 12)
    expect_equal(wide$L2_3[1], 102)
    expect_equal(as.character(wide$event_1[4]), "A")
    expect_equal(wide$L1_1[4], 41)
    expect_equal(wide$L2_1[4], 401)
})

test_that("event tables reject calendar dates", {
    fixture <- event_table_fixture()
    fixture$events$time <- as.Date("2026-01-01") + fixture$events$time
    combined <- event_table_combined(fixture$baseline, fixture$events, time_covariates = c("L1", "L2"))
    separate <- event_table_separate(fixture$baseline, fixture$events)
    expect_error(prepare_rtmle_data(combined, min_events = 0), "numeric.*since baseline")
    expect_error(prepare_rtmle_data(separate, min_events = 0), "numeric.*since baseline")
})

test_that("event tables reject invalid roles, columns and event values", {
    fixture <- event_table_fixture()
    x <- event_table_init(fixture$baseline)
    events <- fixture$events
    expect_error(add_continuous_data(x, event_data = events, covariate_events = c("L", "visit")),
                 "treatment_events")
    expect_error(add_continuous_data(x, event_data = events, treatment_events = c("A", "review")),
                 "covariate_events")
    for (column in c("id", "time", "event", "A", "L1")) {
        bad <- events
        bad[[column]] <- NULL
        expect_error(event_table_combined(fixture$baseline, bad, time_covariates = c("L1", "L2")),
                     info = paste("missing column", column))
    }
    bad <- events
    bad$event[1] <- "unknown"
    expect_error(event_table_combined(fixture$baseline, bad, time_covariates = c("L1", "L2")),
                 "event|unknown")
    bad <- events
    bad$A[1] <- 2
    expect_error({
        y <- event_table_combined(fixture$baseline, bad, time_covariates = c("L1", "L2"))
        prepare_rtmle_data(y, min_events = 0)
    }, "binary|0/1|treatment")
    bad <- events
    bad$L2[2] <- NA
    expect_error({
        y <- event_table_combined(fixture$baseline, bad, time_covariates = c("L1", "L2"))
        prepare_rtmle_data(y, min_events = 0)
    }, "[Mm]issing|NA|[Nn]onfinite")
    bad <- rbind(events, events[events$event == "Y", ])
    expect_error({
        y <- event_table_combined(fixture$baseline, bad, time_covariates = c("L1", "L2"))
        prepare_rtmle_data(y, min_events = 0)
    }, "[Dd]uplicat|tie")
})

# Every person has a treatment opportunity after an L measurement, with the
# previous treatment determined by their own baseline value.
event_table_switch_fixture <- function(n = 20L, switches = n) {
    baseline_A <- rep(c(0, 1), length.out = n)
    baseline <- data.frame(id = seq_len(n), A = baseline_A,
                           L1 = seq_len(n), L2 = seq_len(n) + 100)
    assigned <- baseline_A
    assigned[seq_len(switches)] <- 1 - baseline_A[seq_len(switches)]
    events <- data.frame(
        id = rep(seq_len(n), each = 3),
        time = rep(c(0.5, 1, 2), n),
        event = rep(c("L", "A", "Y"), n),
        A = as.vector(rbind(rep(NA_real_, n), assigned, rep(NA_real_, n))),
        L1 = as.vector(rbind(seq_len(n), rep(NA_real_, n), rep(NA_real_, n))),
        L2 = as.vector(rbind(seq_len(n) + 100, rep(NA_real_, n), rep(NA_real_, n)))
    )
    list(baseline = baseline, events = events)
}

event_table_prepare_either <- function(fixture, separate, ...) {
    x <- if (separate) event_table_separate(fixture$baseline, fixture$events) else
        event_table_combined(fixture$baseline, fixture$events, time_covariates = c("L1", "L2"))
    prepare_rtmle_data(x, min_events = 0, ...)
}

event_table_check_static <- function(fixture, separate, ...) {
    x <- event_table_prepare_either(fixture, separate, ...)
    continuous_check_static_support(x$continuous_treatment_support)
    x
}

test_that("always switching treatment is rejected for both input formats", {
    fixture <- event_table_switch_fixture()
    for (separate in c(FALSE, TRUE)) {
        expect_no_error(event_table_prepare_either(fixture, separate))
        expect_error(event_table_check_static(fixture, separate), "switch|support",
                     info = paste("separate input:", separate))
    }
})

test_that("nearly universal switching warns but unchanged assignments remain allowed", {
    for (separate in c(FALSE, TRUE)) {
        # 19/20 switches is exactly the 95% warning threshold.
        near <- event_table_switch_fixture(switches = 19)
        expect_warning(event_table_check_static(near, separate), "switch|support",
                       info = paste("separate input:", separate))
        # At 18/20 switches there is no warning or error.
        some_unchanged <- event_table_switch_fixture(switches = 18)
        expect_no_warning(event_table_check_static(some_unchanged, separate))
    }
})

test_that("switching is evaluated relative to each person's previous treatment", {
    fixture <- event_table_switch_fixture()
    first_A <- fixture$events[fixture$events$event == "A", ]
    # A second opportunity repeats the first assignment, so only half of all
    # opportunities switch despite everyone switching at their first one.
    first_A$time <- 1.5
    fixture$events <- rbind(fixture$events, first_A)
    fixture$events <- fixture$events[order(fixture$events$id, fixture$events$time), ]
    for (separate in c(FALSE, TRUE))
        expect_no_warning(event_table_check_static(fixture, separate))
})

test_that("no treatment opportunities are allowed without switch warnings", {
    fixture <- event_table_switch_fixture()
    fixture$events <- fixture$events[fixture$events$event != "A", ]
    for (separate in c(FALSE, TRUE))
        expect_no_warning(event_table_check_static(fixture, separate))
})


test_that("static support uses all observed opportunities even after event-index cutoff", {
    fixture <- event_table_switch_fixture()
    for (separate in c(FALSE, TRUE)) {
        # The first treatment opportunity is dropped by this cutoff. It must
        # still contribute to the support diagnostic for the observed data.
        expect_error(event_table_check_static(fixture, separate, last_non_terminal_event = 0),
                     "switch|support", info = paste("separate input:", separate))
    }
})

test_that("a static analysis checks for all-switching treatment opportunities", {
    fixture <- event_table_switch_fixture()
    for (separate in c(FALSE, TRUE)) {
        x <- event_table_prepare_either(fixture, separate)
        x <- regime(x, "Always_A", intervention = 1, treatment_variables = "A")
        x <- target(x, "risk", regimes = "Always_A")
        x <- model_formula(x, verbose = FALSE)
        expect_error(run_rtmle(x, estimator = "one-step"), "switch|support",
                     info = paste("separate input:", separate))
    }
})

test_that("configured terminal event labels are accepted in event tables", {
    fixture <- event_table_fixture()
    fixture$events$event <- ifelse(fixture$events$event == "Y", "stroke",
                            ifelse(fixture$events$event == "D", "death",
                            ifelse(fixture$events$event == "C", "dropout", fixture$events$event)))
    x <- rtmle_init(c(0, 5), name_id = "id", name_time = "time", name_outcome = "stroke",
                    name_competing = "death", name_censoring = "dropout", continuous_time = TRUE)
    x <- add_baseline_data(x, fixture$baseline)
    x <- add_continuous_data(x, event_data = fixture$events,
                             treatment_events = c("A", "review"),
                             covariate_events = c("L", "visit"), time_covariates = c("L1", "L2"))
    wide <- prepare_rtmle_data(x, min_events = 0)$prepared_data
    expect_equal(as.character(wide$event_2[3]), "C")
    expect_equal(as.character(wide$event_3[c(2, 4)]), c("D", "tauend"))
    expect_equal(as.character(wide$event_4[1]), "Y")
})

test_that("terminal labels cannot also identify treatment or measurement events", {
    fixture <- event_table_fixture()
    x <- event_table_init(fixture$baseline)
    expect_error(add_continuous_data(x, event_data = fixture$events,
                                    treatment_events = c("A", "Y"),
                                    covariate_events = c("L", "visit")), "[Tt]erminal")
    expect_error(add_continuous_data(x, event_data = fixture$events,
                                    treatment_events = c("A", "review"),
                                    covariate_events = c("L", "C")), "[Tt]erminal")
})

test_that("event tables can contain no time-varying covariates or measurements", {
    fixture <- event_table_fixture()
    events <- fixture$events[!fixture$events$event %in% c("L", "visit"), c("id", "time", "event", "A")]
    x <- add_continuous_data(event_table_init(fixture$baseline), event_data = events,
                             treatment_events = c("A", "review"),
                             covariate_events = NULL, time_covariates = character())
    wide <- prepare_rtmle_data(x, min_events = 0)$prepared_data
    expect_equal(as.character(wide$event_1), c("A", "A", "C", "A"))
    expect_false(any(grepl("^L[12]_[0-9]+$", names(wide))))
})

test_that("a single call cannot mix event-table and separate data inputs", {
    fixture <- event_table_fixture()
    common <- list(x = event_table_init(fixture$baseline), event_data = fixture$events,
                   treatment_events = c("A", "review"), covariate_events = c("L", "visit"),
                   time_covariates = c("L1", "L2"))
    examples <- list(outcome_data = data.frame(id = 1, time = 4),
                     censored_data = data.frame(id = 3, time = 4),
                     competing_data = data.frame(id = 2, time = 4),
                     timevar_data = list(A = data.frame(id = 1, time = 1, value = 1)),
                     L3 = data.frame(id = 1, time = 2, value = 1))
    for (part in names(examples)) {
        args <- common
        args[[part]] <- examples[[part]]
        expect_error(do.call(add_continuous_data, args), "cannot be combined", info = part)
    }
})

test_that("an all-switching dynamic rule is not blocked by the static support check", {
    fixture <- event_table_switch_fixture()
    # Preserve observed baseline treatment and switch at each later opportunity.
    rule <- function(data, intervention_table, time_node) {
        if (time_node > 0) {
            data[[paste0("A_", time_node)]] <- 1 - data[[paste0("A_", time_node - 1L)]]
        }
        data
    }
    for (separate in c(FALSE, TRUE)) {
        x <- event_table_prepare_either(fixture, separate)
        x <- regime(x, "switch", intervention = 1, treatment_variables = "A",
                    intervene_function = rule)
        x <- target(x, "risk", regimes = "switch")
        x <- model_formula(x, verbose = FALSE)
        result <- suppressWarnings(run_rtmle(x, learner = "learn_oipcw",
                                              propensity_learner = "learn_glm", estimator = "one-step"))
        expect_equal(result$estimate$Main_analysis$Estimate, 1, tolerance = 1e-7)
    }
})

test_that("event tables honor configured subject and time column names", {
    fixture <- event_table_fixture()
    expected <- prepare_rtmle_data(event_table_combined(fixture$baseline, fixture$events,
                                                        time_covariates = c("L1", "L2")), min_events = 0)
    names(fixture$baseline)[names(fixture$baseline) == "id"] <- "person"
    names(fixture$events)[names(fixture$events) == "id"] <- "person"
    names(fixture$events)[names(fixture$events) == "time"] <- "visit_time"
    x <- rtmle_init(c(0, 5), name_id = "person", name_time = "visit_time",
                    name_outcome = "Y", name_competing = "D", name_censoring = "C",
                    continuous_time = TRUE)
    x <- add_baseline_data(x, fixture$baseline)
    x <- add_continuous_data(x, event_data = fixture$events,
                             treatment_events = c("A", "review"),
                             covariate_events = c("L", "visit"), time_covariates = c("L1", "L2"))
    result <- prepare_rtmle_data(x, min_events = 0)
    expect_equal(result$continuous_prepared$wide_data, expected$continuous_prepared$wide_data)
    expect_equal(result$prepared_data$person, 1:4)
})

test_that("simulated event tables agree with the existing separate-input workflow", {
    set.seed(34)
    d <- simulate_continuous_time_data(1000, uncensored = FALSE, no_competing_events = FALSE)
    x <- rtmle_init(c(0, 720), name_id = "id", name_time = "time", name_outcome = "Y",
                    name_competing = "D", name_censoring = "C", continuous_time = TRUE)
    x <- add_baseline_data(x, d$baseline_data)
    combined <- add_continuous_data(x, event_data = d$timevarying_data,
                                    treatment_events = "A", covariate_events = "L", time_covariates = "L")
    separate <- continuous_add_example_events(x, d)
    combined <- prepare_rtmle_data(combined)
    separate <- prepare_rtmle_data(separate)
    expect_equal(combined$continuous_prepared$wide_data, separate$continuous_prepared$wide_data)
    expect_equal(combined$continuous_prepared$info, separate$continuous_prepared$info)
    fit <- function(x) {
        x <- regime(x, "Always_A", intervention = 1, treatment_variables = "A")
        x <- target(x, "risk", regimes = "Always_A")
        x <- model_formula(x, verbose = FALSE)
        run_rtmle(x, learner = "learn_oipcw", propensity_learner = "learn_glm", estimator = "one-step")
    }
    a <- fit(combined)
    b <- fit(separate)
    expect_equal(a$estimate$Main_analysis, b$estimate$Main_analysis, tolerance = 1e-10)
    expect_equal(a$IC$risk$Always_A, b$IC$risk$Always_A, tolerance = 1e-10)
})


test_that("event tables cannot use the treatment column as the identifier or time", {
    for (kind in c("id", "time")) {
        x <- rtmle_init(c(0, 2), name_id = if (kind == "id") "A" else "id",
                        name_time = if (kind == "time") "A" else "time",
                        name_outcome = "Y", continuous_time = TRUE)
        data <- data.frame(id = 1, time = 1, A = 0, event = "Y")
        expect_error(add_continuous_data(x, event_data = data,
                                         treatment_events = "treat", covariate_events = NULL),
                     "distinct")
    }
})


test_that("a replacement event table does not retain old covariate streams", {
    fixture <- event_table_fixture()
    x <- event_table_combined(fixture$baseline, fixture$events, time_covariates = c("L1", "L2"))
    events <- fixture$events[!fixture$events$event %in% c("L", "visit"), c("id", "time", "event", "A")]
    x <- add_continuous_data(x, event_data = events,
                             treatment_events = c("A", "review"), covariate_events = NULL,
                             time_covariates = character())
    expect_equal(names(x$long_data$timevar_data), "A")
    expect_no_error(prepare_rtmle_data(x, min_events = 0))
})
