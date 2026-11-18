test_that("continuous TMLE solves the same weighted intercept score it applies", {
    Y <- c(0, 0.4, 1.5, 0.2)
    Q <- c(0.15, 0.3, 0.6, 0.8)
    weights <- c(1, 2, 1, 1)
    updated <- continuous_tmle_update(Y, Q, weights)

    expect_equal(updated$predicted_values,
                 stats::plogis(stats::qlogis(Q) + updated$epsilon))
    expect_equal(sum(weights * (Y - updated$predicted_values)), 0,
                 tolerance = 1e-8)
    expect_equal(updated$score, 0, tolerance = 1e-8)
    expect_true(all(updated$predicted_values >= 0 &
                    updated$predicted_values <= 1))
})

test_that("continuous TMLE preserves predictions without positive weights", {
    Q <- c(0, 0.3, 1)
    updated <- continuous_tmle_update(c(0, 2, 1), Q, rep(0, 3))
    expect_identical(updated$predicted_values, Q)
    expect_equal(updated$epsilon, 0)
    expect_equal(updated$score, 0)
})

test_that("continuous TMLE uses explicit limits at boundary means", {
    zero <- continuous_tmle_update(c(0, 0, 5), c(0, 0.3, 1), c(1, 2, 0))
    one <- continuous_tmle_update(c(1, 1, 5), c(0, 0.3, 1), c(1, 2, 0))
    expect_equal(zero$predicted_values, rep(0, 3))
    expect_equal(zero$epsilon, -Inf)
    expect_equal(one$predicted_values, rep(1, 3))
    expect_equal(one$epsilon, Inf)
    expect_equal(zero$score, 0)
    expect_equal(one$score, 0)
})

test_that("continuous TMLE tolerates floating-point error at the upper boundary", {
    weights <- c(0.00110936741584863, 1191.8191362854075, 0.0037395730344208732,
                 5.6950401078635444, 0.0126777718939677323, 0.0004569665949320575)
    one <- continuous_tmle_update(rep(1, 6), rep(0.3, 6), weights)
    expect_equal(one$predicted_values, rep(1, 6))
    expect_equal(one$epsilon, Inf)
    expect_equal(one$score, 0)
    rounded <- continuous_tmle_update(rep(1 + .Machine$double.eps, 6),
                                      rep(0.3, 6), weights)
    expect_equal(rounded$predicted_values, rep(1, 6))
    expect_error(continuous_tmle_update(rep(1 + 1e-10, 6), rep(0.3, 6), weights),
                 "no solution")
})

test_that("continuous TMLE rejects invalid inputs and infeasible scores", {
    expect_error(continuous_tmle_update(c(2, 2), c(0.2, 0.3), c(1, 1)),
                 "no solution")
    expect_error(continuous_tmle_update(c(0, -1), c(0.2, 0.3), c(1, 1)),
                 "non-negative")
    expect_error(continuous_tmle_update(c(0, 1), c(0.2, 0.3), c(1, -1)),
                 "non-negative")
    expect_error(continuous_tmle_update(c(0, 1), c(0.2, NA), c(1, 1)),
                 "finite")
    expect_error(continuous_tmle_update(c(0, 1), c(0.2, 1.3), c(1, 1)),
                 "\\[0, 1\\]")
    expect_error(continuous_tmle_update(c(0, 1), 0.2, c(1, 1)),
                 "same length")
    expect_error(continuous_tmle_update(0, 0.2, 1, c(0, 1)),
                 "prediction_range")
})

test_that("the IPCW learner fits responses exceeding one and aliased predictors", {
    d <- data.frame(Y = c(0, 1.2, 0, 1.2, 0.2, 1, 0.2, 1),
                    X = rep(0:1, each = 4))
    d$alias <- 2 * d$X
    prediction_data <- data.frame(X = c(0, 0.3, 1), alias = c(0, 0.6, 2))
    fitted <- learn_oipcw("Y ~ X + alias", d, prediction_data,
                          save_fitted_objects = TRUE)
    expect_equal(fitted$predicted_values, rep(0.6, 3), tolerance = 1e-8)
    reused <- learn_oipcw("Y ~ X + alias", NULL, prediction_data,
                          reuse_fit = fitted)
    expect_equal(reused$predicted_values, fitted$predicted_values)
})

test_that("the IPCW learner agrees with logistic regression on bounded responses", {
    set.seed(91)
    d <- data.frame(X = stats::rnorm(80))
    d$Y <- stats::rbinom(80, 1, stats::plogis(0.2 + 0.4 * d$X))
    expected <- stats::predict(stats::glm(Y ~ X, d, family = stats::binomial()),
                               newdata = d, type = "response")
    fitted <- learn_oipcw("Y ~ X", d, d)
    expect_equal(fitted$predicted_values, unname(expected), tolerance = 1e-7)
})

test_that("the IPCW learner handles constants and rejects invalid responses", {
    d <- data.frame(Y = rep(0, 4), X = factor(rep("only", 4)))
    expect_equal(learn_oipcw("Y ~ X", d, d)$predicted_values, rep(0, 4))
    d$Y <- rep(1, 4)
    expect_equal(learn_oipcw("Y ~ X", d, d)$predicted_values, rep(1, 4))
    d$Y <- c(0, 1.2, 0, 1.2)
    expect_equal(learn_oipcw("Y ~ X", d, d)$predicted_values,
                 rep(0.6, 4), tolerance = 1e-8)
    d$Y <- c(0, -1, 0, 1)
    expect_error(learn_oipcw("Y ~ X", d, d), "non-negative")
    d$Y <- rep(2, 4)
    expect_error(learn_oipcw("Y ~ X", d, d), "no finite solution")
})
