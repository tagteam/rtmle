library(testthat)
library(rtmle)

test_that("treatment trajectories show continuous and discrete representations", {
    data(rtmle_object)

    plot <- plot_treatment_trajectories(rtmle_object)
    expect_s3_class(plot, "ggplot")
    expect_setequal(
        unique(as.character(plot$data$panel)),
        c("Continuous time", "Discrete time")
    )
    expect_equal(max(plot$data$subject), NROW(rtmle_object$prepared_data))
    expect_equal(
        plot$scales$get_scales("y")$limits,
        c(1, NROW(rtmle_object$prepared_data) + 1)
    )
    expect_identical(plot$coordinates$clip, "off")
    expect_true(all(c("A=0, B=0", "A=1, B=0", "A=0, B=1", "A=1, B=1") %in%
                    levels(plot$data$state)))
    expect_true(inherits(plot$theme$axis.text.y, "element_blank"))

    built <- ggplot2::ggplot_build(plot)
    expect_length(built$data, 3L)
    expect_true(any(vapply(
        plot$layers,
        function(layer) inherits(layer$geom, "GeomPoint"),
        logical(1L)
    )))
    expect_gt(nrow(built$data[[3L]]), 0L)
    expect_setequal(
        unique(built$data[[1L]]$xintercept),
        rtmle_object$time_grid_scale
    )
    expect_silent(ggplot2::ggplot_gtable(built))
})

test_that("treatment trajectories use whichever representation is available", {
    data(rtmle_object)

    continuous_object <- rtmle_object
    continuous_object$data <- list()
    continuous_object$prepared_data <- NULL
    continuous_object$followup <- NULL
    continuous_plot <- plot_treatment_trajectories(
        continuous_object,
        treatment_variables = c("A", "B")
    )
    expect_identical(
        unique(as.character(continuous_plot$data$panel)),
        "Continuous time"
    )
    expect_silent(ggplot2::ggplot_build(continuous_plot))

    discrete_object <- rtmle_object
    discrete_object$long_data <- NULL
    discrete_plot <- plot_treatment_trajectories(
        discrete_object,
        treatment_variables = c("A", "B")
    )
    expect_identical(
        unique(as.character(discrete_plot$data$panel)),
        "Discrete time"
    )
    expect_silent(ggplot2::ggplot_build(discrete_plot))
})

test_that("exposure-interval long data treat absent intervals as unexposed", {
    data(simulated_cohort)
    registered <- register_format(
        simulated_cohort,
        treatment_variables = c("A", "B")
    )
    object <- rtmle_init(
        time_grid = seq(0, 20, 4),
        name_id = "id",
        name_outcome = "stroke",
        name_competing = "death",
        name_censoring = "dropout",
        censored_label = "censored"
    )
    object <- add_long_data(
        object,
        outcome_data = registered$timevar_data$stroke[!duplicated(id)],
        censored_data = registered$timevar_data$dropout,
        competing_data = registered$timevar_data$death,
        timevar_data = registered$timevar_data[c("A", "B")]
    )
    object <- add_baseline_data(object, registered$baseline_data)

    plot <- plot_treatment_trajectories(
        object,
        treatment_variables = c("A", "B")
    )
    expect_false(any(grepl("Missing", levels(plot$data$state), fixed = TRUE)))
    expect_silent(ggplot2::ggplot_build(plot))
})

test_that("continuous trajectories show observed points and put missing last", {
    x <- rtmle_init(
        time_grid = c(0, 2, 4),
        name_id = "id",
        name_outcome = "Y",
        name_competing = NULL,
        name_censoring = "Censored",
        censored_label = "censored"
    )
    x <- add_long_data(
        x,
        outcome_data = data.frame(id = 1, date = 4),
        timevar_data = list(
            LDL = data.frame(
                id = c(1, 1, 2),
                date = c(1, 3, 2),
                value = c(120, 130, 110)
            ),
            statins = data.frame(
                id = 1,
                start_date = 0.5,
                end_date = 1.5
            )
        )
    )

    plot <- plot_treatment_trajectories(x, treatment_variables = "LDL")
    built <- ggplot2::ggplot_build(plot)
    point_data <- built$data[[length(built$data)]]

    expect_equal(sort(point_data$x), c(1, 2, 3))
    expect_equal(tail(levels(plot$data$state), 1L), "LDL=Missing")

    statin_plot <- plot_treatment_trajectories(
        x,
        treatment_variables = "statins"
    )
    statin_points <- ggplot2::ggplot_build(statin_plot)$data[[2L]]
    expect_equal(sort(statin_points$x), c(0.5, 1.5))
})
