library(testthat)
library(data.table)
library(rtmle)

coef.plot_model_coefficients_test <- function(object, ...) {
    data.table(
        time = c("time_10", "time_2", "time_10", "time_2"),
        node = "n",
        outcome = c("o10", "o2", "o10b", "o2b"),
        terms = "z",
        beta = 1:4
    )
}
registerS3method(
    "coef",
    "plot_model_coefficients_test",
    coef.plot_model_coefficients_test,
    envir = asNamespace("stats")
)

test_that("model coefficient time nodes are ordered numerically", {
    x <- structure(
        list(
            learner = list(fun = "learn_glm", name = "test"),
            models = list(time_2 = list(), time_10 = list()),
            run_time_horizons = 1,
            time_grid = c(0, 1),
            regimes = list(test_regime = list())
        ),
        class = c("plot_model_coefficients_test", "rtmle")
    )

    output <- capture.output(result <- plot_model_coefficients(
        x,
        time_horizon = 1,
        regime = "test_regime",
        times = c("time_10", "time_2"),
        nodes = "n",
        plot_style = "manhattan",
        show_x_labels = TRUE
    ))
    expect_s3_class(result, "ggplot")
    expect_equal(
        unique(as.character(result$data$x_cat)),
        c(
            "time_2 | n | o2",
            "time_2 | n | o2b",
            "time_10 | n | o10",
            "time_10 | n | o10b"
        )
    )
})
