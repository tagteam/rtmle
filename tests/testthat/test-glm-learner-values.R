test_that("logistic fits recover the two observed group risks", {
    # Each group contains 40 observations, with 10 and 30 events respectively.
    d <- data.frame(Y = rep(c(1, 0, 1, 0), c(10, 30, 30, 10)),
                    A_0 = rep(0:1, each = 40), AA = rep(1:5, 16))
    new <- d[c(41, 1, 41), ]
    fit <- learn_glm("Y ~ A_0 + AA", d, new, learn_variables = "A",
                     save_fitted_objects = TRUE)
    expect_equal(as.numeric(fit$predicted_values), c(.75, .25, .75))
    expect_equal(unname(fit$fit_summary[, 1]),
                 c(qlogis(.25), qlogis(.75) - qlogis(.25)), tolerance = 1e-6)
    reused <- learn_glm("Y ~ A_0", d, new[2:1, ], reuse_fit = fit)
    expect_equal(as.numeric(reused$predicted_values), c(.25, .75))

    d$Y[1] <- NA_real_
    constant <- learn_glm("Y ~ A_0", d, new, learn_variables = "unused",
                          save_fitted_objects = TRUE)
    expect_equal(constant$predicted_values, rep(39/79, 3))
    expect_equal(as.numeric(constant$fit), qlogis(39/79))
    expect_equal(learn_glm("Y ~ A_0", d, new[1, ], reuse_fit = constant)$predicted_values,
                 39/79)
})

test_that("intercept-only logistic fits predict the event fraction and reject invalid outcomes", {
    d <- data.frame(Y = rep(c(0, 0, 0, 1), 12), A = rep(0:1, 24))
    fit <- learn_glm("Y ~ 1", d, data.frame(A = 0:1),
                     save_fitted_objects = TRUE)
    expect_equal(as.numeric(fit$predicted_values), c(.25, .25))
    expect_equal(unname(fit$fit_summary[1, 1]), qlogis(.25), tolerance = 1e-6)
    d$Y[1] <- 2
    expect_error(learn_glm("Y ~ 1", d, d), "Could not fit model with glm")
})

test_that("Gaussian lasso recovers the soft-threshold solution and keeps prediction order", {
    # With centered unit-variance z and slope 3, the lasso slope is 3 - lambda.
    d <- data.frame(z = rep(c(-1, 1), 20))
    d$Y <- 2 + 3 * d$z + rep(c(-.2, -.2, .2, .2), 10)
    d <- rbind(d, data.frame(z = 9, Y = NA_real_))
    new <- data.frame(z = c(1, NA_real_, -1, 0))
    fit <- learn_glmnet("Y ~ z", d, new, lambda = c(1, .5),
                        save_fitted_objects = TRUE)
    expect_equal(fit$fit_summary$selected_lambda, .5)
    expect_equal(as.numeric(fit$fit_summary$selected_beta["z", ]), 2.5)
    expect_equal(fit$predicted_values, c(4.5, NA_real_, -.5, 2))
    reused <- learn_glmnet("Y ~ z", d, new[4:1, , drop = FALSE], reuse_fit = fit)
    expect_equal(reused$predicted_values, rev(fit$predicted_values))
    expect_error(learn_glmnet("Y ~ z", d, new, selector = "largest"), "selector")
})

test_that("binary glmnet and its cross-validation selectors recover known risks", {
    d <- data.frame(Y = rep(c(1, 0, 1, 0), c(10, 30, 30, 10)),
                    z = rep(c(-1, 1), each = 40))
    new <- data.frame(z = c(1, -1))
    fit <- learn_glmnet("Y ~ z", d, new, lambda = 0)
    expect_equal(fit$predicted_values, c(.75, .25), tolerance = 1e-6)
    # Each fold has eight observations from each group. The event counts vary
    # across folds so the minimum and one-standard-error rules select different penalties.
    foldid <- c(rep(1:5, 0:4), rep(1:5, 8:4), rep(1:5, 8:4), rep(1:5, 0:4))
    lambda <- c(.22, .10, .01)
    reference <- glmnet::cv.glmnet(stats::model.matrix(~ z, d), d$Y, family = "binomial",
                                   lambda = lambda, foldid = foldid)
    expect_gt(reference$lambda.1se, reference$lambda.min)
    for (selection in c("min", "1se")) {
        penalty <- reference[[paste0("lambda.", selection)]]
        fit <- learn_glmnet("Y ~ z", d, new, selector = selection,
                            lambda = lambda, foldid = foldid)
        expect_equal(fit$fit_summary$selected_lambda, penalty)
        expect_equal(as.numeric(fit$fit_summary$selected_beta["z", ]),
                     as.numeric(stats::coef(reference, s = penalty)["z", ]))
        # Logistic score equation for standardized z = +/-1 gives .75 - lambda.
        expect_equal(fit$predicted_values, c(.75 - penalty, .25 + penalty), tolerance = 1e-6)
    }
})

test_that("learn_glmnet matches a direct elastic-net fit for binomial predictions", {
    set.seed(921)
    n <- 120L
    d <- data.frame(z = rnorm(n), w = rnorm(n),
                    group = factor(rep(c("a", "b", "c"), length.out = n)))
    d$Y <- rbinom(n, 1, plogis(.2 + d$z - .7 * d$w +
                               .9 * (d$group == "b") - .8 * (d$group == "c")))
    newdata <- data.frame(z = c(1.5, NA, -.8, .2, 1.5),
                          w = c(-.4, .6, 1, -.9, -.4),
                          group = factor(c("a", "c", "b", "c", "a"),
                                         levels = levels(d$group)),
                          unused = NA_real_)
    d$Y[3] <- NA_integer_
    d$w[5] <- NA_real_
    contrasts(d$group) <- stats::contr.treatment(levels(d$group))
    contrasts(newdata$group) <- contrasts(d$group)
    lambda <- .05
    alpha <- .4

    # Build the reference predictors directly from the factor labels.
    design <- function(data) {
        cbind("(Intercept)" = 1, z = data$z, w = data$w,
              groupb = as.numeric(data$group == "b"),
              groupc = as.numeric(data$group == "c"))
    }
    train <- d[-c(3, 5), ]
    reference <- glmnet::glmnet(design(train), train$Y, family = "binomial",
                                lambda = lambda, alpha = alpha)
    rows <- c(1, 3, 4, 5)
    expected <- rep(NA_real_, nrow(newdata))
    expected[rows] <- as.numeric(predict(reference, newx = design(newdata[rows, ]),
                                         s = lambda, type = "response"))
    expect_gt(diff(range(expected, na.rm = TRUE)), .1)

    fit <- learn_glmnet("Y ~ z + w + group", d, newdata,
                        lambda = lambda, alpha = alpha, family = "binomial")
    expect_equal(fit$predicted_values, expected, tolerance = 1e-7)
})

test_that("glm_net maps binary factor levels and validates prediction penalties", {
    d <- data.frame(Y = factor(rep(c("yes", "no", "yes", "no"), c(10, 30, 30, 10)),
                               levels = c("no", "yes")),
                    z = rep(c(-1, 1), each = 40))
    new <- data.frame(z = c(-1, 1))
    fit <- glm_net("Y ~ z", d, lambda = c(.3, .1, .01), cv = FALSE)
    predict_risk <- function(...) riskRegression::predictRisk(fit, newdata = new, ...)
    expect_error(predict_risk(), "single numeric lambda")
    expect_error(predict_risk(lambda = c(.3, .1)), "single numeric lambda")
    expect_error(predict_risk(lambda = "0.1"), "single numeric lambda")
    expect_error(predict_risk(lambda = .02), "not fitted with the specified")
    expect_error(predict_risk(lambda = .1, times = 1:2), "wrong dimensions")
    attr(fit, "selected.lambda") <- .1
    expect_equal(as.numeric(predict_risk()), c(.35, .65), tolerance = 1e-6)
    expect_equal(as.numeric(predict_risk(lambda = .01)), c(.26, .74), tolerance = 1e-6)
})

test_that("Gaussian glm_net predicts at an explicit cross-validated penalty", {
    d <- data.frame(z = rep(c(-1, 1), 20))
    d$Y <- 2 + 3 * d$z + rep(c(-.2, -.2, .2, .2), 10)
    fit <- glm_net(Y ~ z, d, lambda = c(1, .5, .25), foldid = rep(1:5, 8))
    expect_true(fit$cv)
    expect_equal(as.numeric(riskRegression::predictRisk(fit, data.frame(z = c(-1, 1)),
                                                         lambda = .5)),
                 c(-.5, 4.5))
})
