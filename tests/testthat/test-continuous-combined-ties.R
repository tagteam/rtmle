combined_tie_example <- function() {
    x <- rtmle_init(c(0, 3), name_id = "id", name_time = "time",
                    name_outcome = "Y", name_competing = "D",
                    name_censoring = "C", continuous_time = TRUE)
    x <- add_baseline_data(x, data.frame(id = 1:3, A = c(0, 1, 0),
                                         L1 = c(10, 20, 30), L2 = c(100, 200, 300)))
    events <- data.frame(
        id = c(1, 1, 1, 2, 2, 3), time = c(1, 1, 3, 1, 3, 1),
        event = c("A", "visit", "Y", "A", "C", "visit"),
        A = c(1, NA, NA, 1, NA, NA),
        L1 = c(NA, 11, 12, 21, 22, 31),
        L2 = c(NA, 101, 102, 201, 202, 301)
    )
    list(x = x, events = events)
}

add_combined_ties <- function(fixture, events = fixture$events) {
    add_continuous_data(fixture$x, event_data = events,
                        treatment_events = "A", covariate_events = "visit",
                        time_covariates = c("L1", "L2"))
}

test_that("tied treatment and covariate rows form one treatment event", {
    fixture <- combined_tie_example()
    x <- add_combined_ties(fixture)
    wide <- prepare_rtmle_data(x, min_events = 0)$prepared_data
    expect_equal(as.character(wide$event_1), c("A", "A", "L"))
    expect_equal(wide$A_1, c(1, 1, 0))
    expect_equal(wide$L1_1, c(11, 21, 31))
    expect_equal(wide$L2_1, c(101, 201, 301))
    # Values on a single terminal row remain part of that terminal event.
    expect_equal(as.character(wide$event_2[1:2]), c("Y", "C"))
    expect_equal(wide$L1_2[1:2], c(12, 22))
    expect_equal(wide$L2_2[1:2], c(102, 202))
})

test_that("same-time nonterminal records must not disagree", {
    fixture <- combined_tie_example()
    duplicate <- fixture$events[1, ]
    duplicate$A <- 0
    expect_error(add_combined_ties(fixture, rbind(fixture$events, duplicate)),
                 "Conflicting values for A.*subject 1.*time 1")

    # A covariate recorded on both rows of a tied A/visit must agree.
    conflicting <- fixture$events
    conflicting$L1[1] <- 99
    expect_error(add_combined_ties(fixture, conflicting),
                 "Conflicting values for L1.*subject 1.*time 1")

    matching <- fixture$events
    matching$L1[1] <- 11
    matching$L2[1] <- 101
    original <- prepare_rtmle_data(add_combined_ties(fixture), min_events = 0)$prepared_data
    result <- prepare_rtmle_data(add_combined_ties(fixture, matching), min_events = 0)$prepared_data
    expect_equal(result, original)

    # Repeated identical nonterminal records describe the same observation.
    repeated <- rbind(fixture$events, fixture$events[2, ])
    result <- prepare_rtmle_data(add_combined_ties(fixture, repeated), min_events = 0)$prepared_data
    expect_equal(result, original)
})

test_that("a terminal row cannot be tied to another event row", {
    fixture <- combined_tie_example()
    for (mark in c("A", "visit", "D", "Y")) {
        extra <- fixture$events[1, ]
        extra$time <- 3
        extra$event <- mark
        if (mark == "visit") {
            extra$L1 <- 12
            extra$L2 <- 102
        }
        expect_error(add_combined_ties(fixture, rbind(fixture$events, extra)),
                     "terminal event.*subject 1.*time 3")
    }
})
