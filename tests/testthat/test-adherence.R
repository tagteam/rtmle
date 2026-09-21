library(testthat)
library(data.table)

make_adherence_fixture <- function() {
    list(
        followup = data.table(
            id = 1:4,
            last_interval = c(2, 2, 1, 2)
        ),
        prepared_data = data.table(
            outcome_1 = c(0, 0, 0, 0),
            outcome_2 = c(0, 0, 0, 0),
            outcome_3 = c(1, 0, 0, 0),
            death_1 = c(0, 0, 0, 0),
            death_2 = c(0, 0, 0, 0),
            death_3 = c(0, 1, 0, 0),
            dropout_0 = c("uncensored", "uncensored", "uncensored", "uncensored"),
            dropout_1 = c("uncensored", "uncensored", "uncensored", "uncensored"),
            dropout_2 = c("uncensored", "uncensored", "censored", "uncensored"),
            dropout_3 = c("uncensored", "uncensored", "uncensored", "censored")
        ),
        names = list(
            outcome = "outcome",
            competing = "death",
            censoring = "dropout",
            censored_label = "censored"
        ),
        regimes = list(
            P1 = list(intervention_match = matrix(
                c(
                    1, 1, 1,
                    1, 0, 0,
                    0, 1, 1,
                    1, 1, 1
                ),
                nrow = 4,
                byrow = TRUE
            )),
            P2 = list(intervention_match = matrix(
                c(
                    1, 0, 0,
                    0, 1, 1,
                    1, 1, 1,
                    1, 1, 1
                ),
                nrow = 4,
                byrow = TRUE
            ))
        ),
        time_grid = 0:2,
        time_grid_labels = as.character(0:2)
    )
}


test_that("adherence collects non-adherence data for selected regimes", {
    x <- make_adherence_fixture()
    non_adherence <- adherence(x)

    expect_s3_class(non_adherence, "data.table")
    expect_equal(
        as.character(non_adherence$regime),
        c("P1", "P1", "P1", "P2", "P2", "P2")
    )
    expect_equal(
        non_adherence$time_nonadherence,
        c(2, 2, 2, 2, 1, 2)
    )
    expect_equal(
        non_adherence$event_nonadherence,
        c(2, 1, 0, 1, 0, 0)
    )
    expect_equal(non_adherence$id, c(1, 2, 4, 1, 3, 4))

    selected <- adherence(x, regimes = "P2")
    expect_equal(as.character(selected$regime), rep("P2", 3L))
    expect_equal(selected$time_nonadherence, c(2, 1, 2))
})


test_that("summary_adherence counts adherence and end-of-follow-up events", {
    x <- make_adherence_fixture()
    adherence_summary <- summary_adherence(x)

    expect_s3_class(adherence_summary, "data.table")
    expect_equal(
        as.character(adherence_summary$regime),
        c("P1", "P1", "P1", "P2", "P2", "P2")
    )
    expect_equal(adherence_summary$time_node, c(0L, 1L, 2L, 0L, 1L, 2L))
    expect_equal(adherence_summary$n_initiated, rep(3L, 6L))
    expect_equal(adherence_summary$n_adherent, rep(c(3L, 2L, 2L), 2L))
    expect_equal(adherence_summary$n_outcome, c(0L, 0L, 1L, 0L, 0L, 1L))
    expect_equal(adherence_summary$n_death, c(0L, 0L, 1L, 0L, 0L, 0L))
    expect_equal(adherence_summary$n_censored, c(0L, 0L, 1L, 0L, 1L, 1L))

    selected <- summary_adherence(x, regimes = "P2")
    expect_equal(as.character(selected$regime), rep("P2", 3L))
    expect_equal(selected$n_censored, c(0L, 1L, 1L))
})


test_that("summary.adherence returns the cached adherence summary", {
    x <- make_adherence_fixture()
    a <- adherence(x, regimes = "P1")

    expect_s3_class(a, "adherence")
    expect_equal(summary(a), summary_adherence(x, regimes = "P1"))
})


test_that("adherence handles objects without censoring variables", {
    x <- make_adherence_fixture()
    x$names$censoring <- NULL
    x$prepared_data <- NULL

    non_adherence <- adherence(x, regimes = "P1")
    expect_equal(non_adherence$time_nonadherence, c(2, 2, 2))
    expect_equal(non_adherence$event_nonadherence, c(2, 1, 2))
})


test_that("plot_adherence uses the selected adherence data", {
    x <- make_adherence_fixture()
    plot <- plot_adherence(x, regimes = "P1")

    expect_s3_class(plot, "ggplot")
})
