library(testthat)

test_that("register_format can return treatment exposure intervals", {
    data(simulated_cohort, package = "rtmle")

    registered <- register_format(
        simulated_cohort,
        treatment_variables = "A"
    )

    expect_identical(
        names(registered$timevar_data$A),
        c("id", "start_date", "end_date")
    )
    expect_identical(
        names(registered$timevar_data$B),
        c("id", "date", "value")
    )
    expect_true(all(registered$timevar_data$A$start_date <=
                    registered$timevar_data$A$end_date))

    ## The first zero-valued observation after an exposed run closes the
    ## interval; a subject with no later treatment observation is closed at
    ## the final observed cohort time.
    history <- simulated_cohort[
        simulated_cohort$id == 5 &
            simulated_cohort$event %in% c("baseline", "visit"),
        c("time", "A")
    ]
    expected_end <- history$time[which(history$A == 0)[1L]]
    observed_end <- registered$timevar_data$A$end_date[
        registered$timevar_data$A$id == 5
    ]
    expect_equal(observed_end[[1L]], expected_end)
})


test_that("register_format keeps the default measurement format", {
    data(simulated_cohort, package = "rtmle")

    registered <- register_format(simulated_cohort)

    expect_identical(
        names(registered$timevar_data$A),
        c("id", "date", "value")
    )
})

test_that("protocol remains available with regime as its alias", {
    x <- rtmle_init(time_grid = 0:1,
                    name_id = "id",
                    name_outcome = "Y")
    intervention <- data.frame(
        time_node = x$intervention_nodes,
        A = factor("1", levels = c("0", "1"))
    )

    x_protocol <- protocol(x, name = "Always_A", intervention = intervention)
    x_regime <- regime(x, name = "Always_A", intervention = intervention)

    expect_identical(x_protocol$regimes, x_regime$regimes)
    expect_null(x_protocol$protocols)
    expect_true(identical(protocol, regime))
})
