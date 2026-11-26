test_that("participant selection rejects ambiguous or empty subsets", {
    data("rtmle_object", package = "rtmle", envir = environment())
    x <- rtmle_object
    unprepared <- rtmle_init(time_grid = 0:1, name_id = "id", name_outcome = "Y")
    expect_error(subset(unprepared, TRUE), "prepare_rtmle_data")
    expect_error(subset(x), "Supply either")
    expect_error(subset(x, TRUE, ids = 1L), "Supply either")
    expect_error(subset(x, 1L), "logical vector")
    expect_error(subset(x, c(TRUE, FALSE)), "one value per participant")
    expect_error(subset(x, NA), "No participants were selected")
    expect_error(subset(x, ids = integer()), "No participants were selected")

    first_id <- x$prepared_data$id[[1L]]
    x$followup <- x$followup[x$followup$id != first_id, ]
    expect_error(subset(x, ids = first_id), "missing follow-up data")
})

test_that("stratification identifies the variable and reports which fit failed", {
    data("rtmle_object", package = "rtmle", envir = environment())
    x <- rtmle_object
    expect_error(stratify(x, absent_column), "variable must be a column")
    x$prepared_data$sex <- as.character(x$prepared_data$sex)
    expect_error(stratify(x, sex), "factor or a numeric vector")
    expect_error(stratify(rtmle_object, sex, targets = "absent target"),
                 "Could not fit stratum.*0.*absent target")
})

test_that("cheap bootstrap validates the fitted estimates and sample sizes", {
    data("rtmle_object", package = "rtmle", envir = environment())
    x <- rtmle_object
    unprepared <- rtmle_init(time_grid = 0:1, name_id = "id", name_outcome = "Y")
    expect_error(cheap_bootstrap(unprepared), "Run run_rtmle")
    expect_error(cheap_bootstrap(x, analyses = "Main_analysis"), "analysis of interest")
    expect_error(cheap_bootstrap(x, time_horizon = 999), "No fitted estimates")
    expect_error(cheap_bootstrap(x, targets = "absent"), "target.*absent")
    expect_error(cheap_bootstrap(x, B = 0), "B > 0")
    expect_error(cheap_bootstrap(x, B = 1.5), "floor")
    expect_error(cheap_bootstrap(x, M = 0), "M > 0")
    expect_error(cheap_bootstrap(x, M = nrow(x$prepared_data)), "smaller than")
    expect_error(cheap_bootstrap(x, B = 3, seeds = 1:2), "length.*seeds")
    x$estimate$Main_analysis$Estimator <- as.character(x$estimate$Main_analysis$Estimator)
    x$estimate$Main_analysis$Estimator[1L] <- "g-formula"
    expect_error(cheap_bootstrap(x), "one estimator")
})

test_that("subsetting a wide-data fit preserves the participant order", {
    data("rtmle_object", package = "rtmle", envir = environment())
    x <- rtmle_object
    # Wide-data analyses have no original long-format records to subset.
    x$long_data <- NULL
    ids <- rev(head(x$prepared_data$id, 4L))
    selected <- subset(x, ids = c(ids, ids[1L]))
    expect_identical(selected$prepared_data$id, rev(ids))
    expect_identical(selected$followup$id, rev(ids))
    expect_null(selected$long_data)
})
