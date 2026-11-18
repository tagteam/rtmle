test_that("GLM learners preserve training factor levels during prediction and reuse", {
    d <- data.table::data.table(
        Y = rep(c(0.1, 0.5, 0.8, 0.4), 10),
        A = rep(c(0, 1, 0, 1), 10),
        event = factor(rep(c("A", "A", "L", "L"), 10),
                       levels = c("A", "L", "C", "Y", "D", "tauend"))
    )
    d[, event_alias := event]
    fitted <- suppressWarnings(learn_glm("Y ~ A + event + event_alias", d, d,
                                         save_fitted_objects = TRUE))
    expect_true(all(is.finite(fitted$predicted_values)))
    reference <- stats::glm(Y ~ A + event, data = d, family = stats::quasibinomial())
    expect_equal(as.numeric(fitted$predicted_values), as.numeric(stats::predict(reference, d, type = "response")),
                 tolerance = 1e-6)
    reused <- suppressWarnings(learn_glm("Y ~ A + event + event_alias", d, d,
                                         reuse_fit = fitted))
    expect_equal(reused$predicted_values, fitted$predicted_values)
    expect_equal(fitted$fit$xlevels$event, c("A", "L"))
})

test_that("the GLM fallback also preserves factor levels", {
    skip_if_not_installed("speedglm")
    testthat::local_mocked_bindings(
        speedglm.wfit = function(...) stop("exercise the fallback"),
        .package = "speedglm"
    )
    d <- data.table::data.table(
        Y = rep(c(0.1, 0.5, 0.8, 0.4), 10),
        event = factor(rep(c("A", "A", "L", "L"), 10),
                       levels = c("A", "L", "C", "Y"))
    )
    fitted <- learn_glm("Y ~ event", d, d, save_fitted_objects = TRUE)
    expect_s3_class(fitted$fit, "glm")
    expect_true(all(is.finite(fitted$predicted_values)))
    prediction <- droplevels(d[event == "A"])
    reused <- learn_glm("Y ~ event", d, prediction, reuse_fit = fitted)
    expect_equal(as.numeric(reused$predicted_values), rep(0.3, nrow(prediction)), tolerance = 1e-6)
})
