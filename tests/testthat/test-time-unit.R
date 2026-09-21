library(testthat)
library(rtmle)

test_that("rtmle_init stores an optional time-unit label", {
    x <- rtmle_init(time_grid = 0:1,
                    name_id = "id",
                    name_outcome = "Y")
    expect_null(x$time_unit)

    x_months <- rtmle_init(time_grid = 0:1,
                           name_id = "id",
                           name_outcome = "Y",
                           time_unit = "Months since treatment start")
    expect_identical(x_months$time_unit, "Months since treatment start")

    expect_error(
        rtmle_init(time_grid = 0:1,
                   name_id = "id",
                   name_outcome = "Y",
                   time_unit = 1),
        "time_unit must be NULL or a non-empty character string"
    )
})

test_that("time-based plots use the configured time-unit label", {
    data(rtmle_object)
    x <- rtmle_object
    x$time_unit <- "Months since treatment start"

    risk_plot <- ggplot2::autoplot(x, conf_int = FALSE)
    adherence_plot <- plot_adherence(x)
    ipw_plot <- plot_IPW(x)

    expect_equal(risk_plot$labels$x, x$time_unit)
    expect_equal(adherence_plot$labels$x, x$time_unit)
    expect_equal(ipw_plot$labels$x, x$time_unit)
})

test_that("time-based plots retain the default label", {
    data(rtmle_object)

    expect_equal(ggplot2::autoplot(rtmle_object, conf_int = FALSE)$labels$x,
                 "Time")
    expect_equal(plot_adherence(rtmle_object)$labels$x, "Time")
    expect_equal(plot_IPW(rtmle_object)$labels$x, "Time")
})

test_that("risk plots can select targets and show bootstrap bands", {
    data(rtmle_object)
    x <- rtmle_object
    main <- data.table::copy(x$estimate$Main_analysis)
    target_labels <- unique(as.character(main$Target))
    bootstrap <- main[, .(
        B = 2L,
        Target,
        Regime,
        Time_horizon,
        Target_parameter,
        Estimator,
        Bootstrap_estimate = Estimate,
        Bootstrap_standard_error = Standard_error,
        Bootstrap_lower = pmax(0, Lower - 0.01),
        Bootstrap_upper = pmin(1, Upper + 0.01)
    )]
    x$estimate$Cheap_bootstrap <- list(Main_analysis = bootstrap)

    risk_plot <- ggplot2::autoplot(
        x,
        targets = target_labels[[1L]],
        conf_int = FALSE,
        bootstrap_conf_int = TRUE
    )

    expect_s3_class(risk_plot, "ggplot")
    expect_equal(unique(as.character(risk_plot$data$Target)), target_labels[[1L]])
    expect_true(any(vapply(risk_plot$layers, function(layer) {
        is.data.frame(layer$data) &&
            all(c("Bootstrap_lower", "Bootstrap_upper") %in% names(layer$data))
    }, logical(1L))))
    expect_silent(ggplot2::ggplot_build(risk_plot))
    expect_error(
        ggplot2::autoplot(
            rtmle_object,
            bootstrap_conf_int = TRUE
        ),
        "no Cheap_bootstrap results"
    )
})

test_that("plot.rtmle returns a composable ggplot object", {
    data(rtmle_object)

    risk_plot <- plot(rtmle_object, conf_int = FALSE) +
        ggplot2::theme(legend.position = "bottom")

    expect_s3_class(risk_plot, "ggplot")
    expect_silent(ggplot2::ggplot_build(risk_plot))
})
