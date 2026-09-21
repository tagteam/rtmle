library(testthat)
library(data.table)

make_treat_unless_data <- function() {
    data.table(
        id = 1:4,
        bleeding_0 = c(0, 1, 0, 1),
        bleeding_1 = c(1, 0, 1, 1),
        surgery_0 = c(0, 0, 1, 0),
        A_0 = factor(c(0, 1, 0, 1), levels = c(0, 1)),
        A_1 = factor(c(0, 0, 1, 1), levels = c(0, 1)),
        B_0 = factor(c(0, 0, 1, 1), levels = c(0, 1)),
        B_1 = factor(c(1, 0, 1, 0), levels = c(0, 1))
    )
}

make_treat_unless_table <- function() {
    data.table(
        time_node = c(0, 1, 0, 1),
        variable = c("A_0", "A_1", "B_0", "B_1"),
        value = factor(c(1, 1, 0, 0), levels = c(0, 1))
    )
}


test_that("treat_unless changes static actions after a contraindication", {
    data <- make_treat_unless_data()
    intervention_table <- make_treat_unless_table()

    intervened <- treat_unless(
        data,
        intervention_table,
        time_node = 1,
        contra_indication = c("bleeding", "surgery"),
        action = c(A = 0, B = 1)
    )

    expect_equal(as.character(intervened$A_0), c("1", "0", "0", "0"))
    expect_equal(as.character(intervened$A_1), c("0", "0", "0", "0"))
    expect_equal(as.character(intervened$B_0), c("0", "1", "1", "1"))
    expect_equal(as.character(intervened$B_1), c("1", "1", "1", "1"))

    ## The default infinite window includes the current node.
    expect_equal(as.character(intervened$A_1[[1L]]), "0")
})


test_that("NA action retains observed treatment after a contraindication", {
    intervened <- treat_unless(
        make_treat_unless_data(),
        make_treat_unless_table(),
        time_node = 1,
        contra_indication = "bleeding",
        action = NA
    )

    expect_equal(as.character(intervened$A_1), c("0", "0", "1", "1"))
    expect_equal(as.character(intervened$B_0), c("0", "0", "0", "1"))
})


test_that("treat_unless can be registered without a wrapper", {
    data <- make_treat_unless_data()
    x <- rtmle_init(
        time_grid = 0:2,
        name_id = "id",
        name_outcome = "Y",
        name_competing = NULL,
        name_censoring = NULL,
        minority_threshold = 0
    )
    x$prepared_data <- data
    x <- regime(
        x,
        name = "A_until_bleeding",
        intervention = data.frame(
            time_node = 0:1,
            A = factor("1", levels = c("0", "1"))
        ),
        intervene_function = treat_unless(
            contra_indication = "bleeding",
            action = 0
        ),
        verbose = FALSE
    )

    evaluated <- rtmle:::evaluate_intervention(
        x$regimes$A_until_bleeding,
        data,
        x$regimes$A_until_bleeding$intervention_table,
        time_node = 1
    )$data
    expect_equal(as.character(evaluated$A_1), c("0", "0", "0", "0"))
})


test_that("numeric and character actions match factor treatment options", {
    for (action in list(0, "0")) {
        data <- make_treat_unless_data()
        data[, c("A_0", "A_1") := lapply(
            .SD,
            as.integer
        ), .SDcols = c("A_0", "A_1")]
        x <- rtmle_init(
            time_grid = 0:2,
            name_id = "id",
            name_outcome = "Y",
            name_competing = NULL,
            name_censoring = NULL,
            minority_threshold = 0
        )
        x$prepared_data <- data
        x <- regime(
            x,
            name = "A_until_bleeding",
            intervention = data.frame(
                time_node = 0:1,
                A = factor("1", levels = c("0", "1"))
            ),
            intervene_function = treat_unless(
                contra_indication = "bleeding",
                action = action
            ),
            verbose = FALSE
        )

        evaluated <- rtmle:::evaluate_intervention(
            x$regimes$A_until_bleeding,
            data,
            x$regimes$A_until_bleeding$intervention_table,
            time_node = 1
        )$data
        expect_equal(evaluated$A_1, c(0L, 0L, 0L, 0L))
    }
})


test_that("lookback_window limits the current-node history window", {
    data <- make_treat_unless_data()
    intervention_table <- make_treat_unless_table()

    current_only <- treat_unless(
        data,
        intervention_table,
        time_node = 1,
        contra_indication = "bleeding",
        action = 0,
        lookback_window = 0
    )
    one_previous_node <- treat_unless(
        data,
        intervention_table,
        time_node = 1,
        contra_indication = "bleeding",
        action = 0,
        lookback_window = 1
    )

    expect_equal(as.character(current_only$A_1), c("0", "1", "0", "0"))
    expect_equal(as.character(one_previous_node$A_1), c("0", "0", "0", "0"))
})


test_that("lookback_window is validated", {
    expect_error(
        treat_unless(
            make_treat_unless_data(),
            make_treat_unless_table(),
            time_node = 1,
            contra_indication = "bleeding",
            action = 0,
            lookback_window = -1
        ),
        "lookback_window"
    )
    expect_error(
        treat_unless(
            make_treat_unless_data(),
            make_treat_unless_table(),
            time_node = 1,
            contra_indication = "bleeding",
            action = 0,
            lookback_window = 1.5
        ),
        "lookback_window"
    )
})
