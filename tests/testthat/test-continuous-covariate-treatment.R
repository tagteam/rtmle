# A covariate event may report the current treatment, but cannot assign it.
covariate_treatment_fixture <- function(baseline_first = TRUE) {
    x <- rtmle_init(c(0, 5), name_id = "id", name_time = "time",
                    name_outcome = "Y", name_censoring = "C", continuous_time = TRUE)
    baseline <- data.frame(id = 1:2, A = c(0, 1), L = c(9, 19))
    events <- data.frame(id = rep(1:2, each = 3), time = rep(c(1, 2, 4), 2),
                         event = rep(c("A", "visit", "Y"), 2),
                         A = c(1, 1, 1, 0, 0, 0),
                         L = c(10, 11, 11, 20, 21, 21))
    if (baseline_first) x <- add_baseline_data(x, data.table::copy(baseline))
    list(x = x, baseline = baseline, events = events)
}

add_covariate_treatment_data <- function(fixture, events = fixture$events) {
    add_continuous_data(fixture$x, event_data = events, treatment_events = "A",
                        covariate_events = "visit", time_covariates = "L")
}

prepare_covariate_treatment_data <- function(x) {
    prepare_rtmle_data(x, min_events = 0)
}

test_that("A observed at covariate events must match the treatment history", {
    fixture <- covariate_treatment_fixture()
    x <- prepare_covariate_treatment_data(add_covariate_treatment_data(fixture))
    expect_equal(x$prepared_data$A_1, c(1, 0))
    expect_equal(x$prepared_data$A_2, c(1, 0))

    changed <- fixture$events
    changed$A[2] <- 0
    x <- add_covariate_treatment_data(fixture, changed)
    expect_error(prepare_covariate_treatment_data(x),
                 "Treatment A cannot change at a covariate event.*subject 1.*time 2")

    changed$A[2] <- NA
    expect_no_error(prepare_covariate_treatment_data(add_covariate_treatment_data(fixture, changed)))
    changed$A[2] <- Inf
    expect_error(add_covariate_treatment_data(fixture, changed), "binary|finite")
})

test_that("baseline can be added after event data and replacement clears old observations", {
    fixture <- covariate_treatment_fixture(baseline_first = FALSE)
    changed <- fixture$events
    changed$A[2] <- 0
    x <- add_covariate_treatment_data(fixture, changed)
    x <- add_baseline_data(x, data.table::copy(fixture$baseline))
    expect_error(prepare_covariate_treatment_data(x), "cannot change at a covariate event")

    # Replacing the event table discards observations from the old table.
    x <- add_continuous_data(x, event_data = fixture$events, treatment_events = "A",
                             covariate_events = "visit", time_covariates = "L")
    expect_no_error(prepare_covariate_treatment_data(x))
})

test_that("replacement covariate-event roles discard stale event observations", {
    fixture <- covariate_treatment_fixture()
    changed <- fixture$events
    changed$A[2] <- 0
    x <- add_covariate_treatment_data(fixture, changed)
    expect_error(prepare_covariate_treatment_data(x), "cannot change at a covariate event")

    # Subject 1's old L event is removed along with its L observation. Subject
    # 2's event remains, so its original A observation is still validated.
    replacement <- data.frame(id = 2, time = 2, value = 21)
    for (role in list("L", data.frame(id = 2, time = 2))) {
        revised <- add_continuous_data(x, timevar_data = list(L = replacement),
                                        covariate_events = role)
        expect_no_error(prepare_covariate_treatment_data(revised))
    }

    # Retaining a conflicting event in the new roles must not disable validation.
    expect_error(prepare_covariate_treatment_data(add_continuous_data(
        x, covariate_events = data.frame(id = 1:2, time = 2))),
        "cannot change at a covariate event")
})

test_that("simultaneous treatment and covariate rows use the treatment assignment", {
    fixture <- covariate_treatment_fixture()
    at_treatment <- fixture$events[1, ]
    at_treatment$event <- "visit"
    events <- rbind(fixture$events, at_treatment)
    x <- prepare_covariate_treatment_data(add_covariate_treatment_data(fixture, events))
    expect_equal(as.character(x$prepared_data$event_1[1]), "A")
    expect_equal(x$prepared_data$A_1[1], 1)

    at_treatment$A <- 0
    events <- rbind(fixture$events, at_treatment)
    expect_error(prepare_covariate_treatment_data(add_covariate_treatment_data(fixture, events)),
                 "Treatment A cannot change at a covariate event.*subject 1.*time 1")
})

test_that("post-terminal covariate records do not silently bypass treatment validation", {
    fixture <- covariate_treatment_fixture()
    later <- fixture$events[2, ]
    later$time <- 5
    later$A <- 0
    events <- rbind(fixture$events, later)
    expect_error(prepare_covariate_treatment_data(add_covariate_treatment_data(fixture, events)),
                 "Treatment A cannot change at a covariate event.*subject 1.*time 5")
    later$A <- 1
    events <- rbind(fixture$events, later)
    expect_no_error(prepare_covariate_treatment_data(add_covariate_treatment_data(fixture, events)))
})

test_that("time-zero covariate treatment is checked against initialized baseline", {
    fixture <- covariate_treatment_fixture()
    at_zero <- fixture$events[2, ]
    at_zero$time <- 0
    at_zero$A <- 0
    at_zero$L <- 9
    events <- rbind(fixture$events, at_zero)
    expect_no_error(prepare_covariate_treatment_data(add_covariate_treatment_data(fixture, events)))

    at_zero$A <- 1
    events <- rbind(fixture$events, at_zero)
    expect_error(prepare_covariate_treatment_data(add_covariate_treatment_data(fixture, events)),
                 "Treatment A cannot change at a covariate event.*subject 1.*time 0")

    # Baseline A still takes precedence over a conflicting time-zero A record.
    zero_treatment <- fixture$events[1, ]
    zero_treatment$time <- 0
    zero_treatment$A <- 1
    zero_treatment$L <- 9
    at_zero$A <- 0
    events <- rbind(fixture$events, zero_treatment, at_zero)
    expect_no_error(prepare_covariate_treatment_data(add_covariate_treatment_data(fixture, events)))
})

test_that("known pre-baseline treatment cannot change at a covariate event", {
    fixture <- covariate_treatment_fixture()
    prior_A <- fixture$events[1, ]
    prior_L <- fixture$events[2, ]
    origin <- 0
    prior_A$time <- origin - 2
    prior_A$A <- 1
    prior_L$time <- origin - 1
    prior_L$A <- 0
    events <- rbind(fixture$events, prior_A, prior_L)
    x <- add_covariate_treatment_data(fixture, events)
    expect_error(prepare_covariate_treatment_data(x),
                 "Treatment A cannot change at a covariate event.*subject 1.*time -1")
    named_role <- add_continuous_data(x, covariate_events = "L")
    expect_error(prepare_covariate_treatment_data(named_role),
                 "Treatment A cannot change at a covariate event.*subject 1.*time -1")

    # A recorded before baseline is compared to the preceding treatment,
    # not to A_0, which is explicitly zero for this subject.
    prior_L$A <- 1
    events <- rbind(fixture$events, prior_A, prior_L)
    expect_no_error(prepare_covariate_treatment_data(
        add_covariate_treatment_data(fixture, events)))

    # There is no known treatment to compare before the first A record.
    prior_L$time <- origin - 3
    prior_L$A <- 0
    events <- rbind(fixture$events, prior_A, prior_L)
    expect_no_error(prepare_covariate_treatment_data(
        add_covariate_treatment_data(fixture, events)))

    # A covariate row itself can establish the first observed state.
    # A later covariate row cannot change it without a treatment event.
    earlier_L <- prior_L
    earlier_L$time <- origin - 4
    earlier_L$A <- 1
    prior_L$time <- origin - 3
    prior_L$A <- 0
    events <- rbind(fixture$events, earlier_L, prior_L)
    expect_error(prepare_covariate_treatment_data(
        add_covariate_treatment_data(fixture, events)),
        "Treatment A cannot change at a covariate event.*subject 1.*time -3")
    prior_L$A <- 1
    events <- rbind(fixture$events, earlier_L, prior_L)
    expect_no_error(prepare_covariate_treatment_data(
        add_covariate_treatment_data(fixture, events)))

    # At a tied time, the treatment assignment is used for the L row.
    prior_A$L <- prior_L$L
    prior_L$time <- origin - 2
    events <- rbind(fixture$events, prior_A, prior_L)
    expect_no_error(prepare_covariate_treatment_data(
        add_covariate_treatment_data(fixture, events)))
    prior_L$A <- 0
    events <- rbind(fixture$events, prior_A, prior_L)
    expect_error(prepare_covariate_treatment_data(
        add_covariate_treatment_data(fixture, events)),
        "Treatment A cannot change at a covariate event.*subject 1.*time -2")

    # A replaced role that removes the old event must not retain its check.
    prior_L$time <- origin - 1
    events <- rbind(fixture$events, prior_A, prior_L)
    x <- add_covariate_treatment_data(fixture, events)
    visits <- fixture$events[fixture$events$event == "visit", c("id", "time")]
    x <- add_continuous_data(x, covariate_events = visits)
    expect_no_error(prepare_covariate_treatment_data(x))
})
