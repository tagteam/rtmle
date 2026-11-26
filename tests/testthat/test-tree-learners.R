test_that("ranger probability forests give the observed proportions in each leaf", {
    skip_if_not_installed("ranger")
    d <- data.table::data.table(x = rep(0:1, each = 4),
                                y = factor(c(0, 0, 0, 1, 0, 1, 1, 1)))
    newdata <- data.table::data.table(x = c(1, NA, 0), unused = NA_real_)
    # Every tree sees the full training sample and splits on its only predictor.
    result <- learn_ranger("y ~ x", d, newdata, save_fitted_objects = TRUE,
                           replace = FALSE, sample.fraction = 1, num.trees = 2,
                           min.node.size = 1, num.threads = 1, seed = 17)
    expect_equal(result$predicted_values, c(3/4, NA, 1/4))
    reused <- learn_ranger("y ~ x", d, newdata[3:1], reuse_fit = result)
    expect_equal(reused$predicted_values, c(1/4, NA, 3/4))
    expect_null(reused[["fit"]])
})

test_that("ranger regression forests give the mean response in each leaf", {
    skip_if_not_installed("ranger")
    d <- data.table::data.table(x = rep(0:1, each = 4),
                                y = c(0, .25, .25, .5, .5, .75, .75, 1))
    result <- learn_ranger("y ~ x", d, d, save_fitted_objects = TRUE,
                           replace = FALSE, sample.fraction = 1, num.trees = 2,
                           min.node.size = 1, num.threads = 1, seed = 17)
    expect_equal(result$predicted_values, rep(c(1/4, 3/4), each = 4))
    pred <- learn_ranger("y ~ x", d, d[c(8, 1)], reuse_fit = result)
    expect_equal(pred$predicted_values, c(3/4, 1/4))
})

test_that("learn_ranger predictions match native probability and regression forests", {
    skip_if_not_installed("ranger")
    iris <- datasets::iris
    training <- data.frame(length = iris$Sepal.Length, species = iris$Species)
    newdata <- data.table::data.table(
        length = c(6.8, NA, 5.1, 6.0, 5.8),
        species = factor(c("virginica", "setosa", "setosa", NA, "versicolor"),
                         levels = levels(iris$Species)),
        unused = NA_real_
    )
    complete_rows <- c(1L, 3L, 5L)
    for (probability in c(TRUE, FALSE)) {
        training$y <- if (probability) factor(as.integer(iris$Sepal.Width > 3)) else iris$Sepal.Width / 5
        reference <- ranger::ranger(y ~ length + species, data = training,
                                    probability = probability, num.trees = 31L,
                                    mtry = 1L, min.node.size = 3L, num.threads = 1L, seed = 17L)
        native <- stats::predict(reference, data = newdata[complete_rows, .(length, species)],
                                 num.threads = 1L)$predictions
        expected <- rep(NA_real_, nrow(newdata))
        expected[complete_rows] <- if (probability) native[, "1"] else native

        fit <- learn_ranger("y ~ length + species", data.table::as.data.table(training), newdata,
                            num.trees = 31L, mtry = 1L, min.node.size = 3L,
                            num.threads = 1L, seed = 17L, save_fitted_objects = TRUE)
        expect_equal(fit$predicted_values, expected)
        row_order <- c(5L, 4L, 1L, 3L, 2L)
        reused <- learn_ranger("y ~ length + species", data.table::as.data.table(training),
                               newdata[row_order], reuse_fit = fit)
        expect_equal(reused$predicted_values, expected[row_order])
    }
})

test_that("one unregularized xgboost tree fits the two observed group means", {
    skip_if_not_installed("xgboost", minimum_version = "3.1.2")
    d <- data.table::data.table(x = rep(0:1, each = 4),
                                y = c(0, 0, 0, 1, 0, 1, 1, 1))
    newdata <- data.table::data.table(x = c(0, NA, 1), unused = NA_real_)
    result <- learn_xgboost("y ~ x", d, newdata, save_fitted_objects = TRUE,
                            xgboost.nrounds = 1L, xgboost.learning_rate = 1,
                            xgboost.reg_lambda = 0, xgboost.max_depth = 1L,
                            xgboost.nthreads = 1L)
    expect_equal(result$predicted_values, c(1/4, NA, 3/4), tolerance = 1e-6)
    reused <- learn_xgboost("y ~ x", d, newdata[3:1], reuse_fit = result)
    expect_equal(reused$predicted_values, c(3/4, NA, 1/4), tolerance = 1e-6)
    expect_null(reused[["fit"]])
})

test_that("learn_xgboost predictions match a direct xgboost fit", {
    skip_if_not_installed("xgboost", minimum_version = "3.1.2")

    train <- data.table::as.data.table(expand.grid(
        x = seq(-2, 2, length.out = 17), group = factor(c("a", "b", "c"))
    ))
    train[, y := plogis(-1 + x + 1.5 * (group == "b") - 2 * x * (group == "c"))]
    newdata <- data.table::data.table(
        group = factor(c("c", "a", "b", "a", "b", "a"), levels = levels(train$group)),
        x = c(1.8, .9, NA, -1.9, .6, .9),
        unused = NA_real_
    )
    train$x[2] <- NA_real_
    train$y[19] <- NA_real_
    seed <- 23L

    # Build the upstream design matrix independently of the wrapper.
    design <- function(d) stats::model.matrix(~ x * group, d)[, -1L, drop = FALSE]
    complete_train <- train[stats::complete.cases(train)]
    set.seed(seed)
    reference <- xgboost::xgboost(
        x = design(complete_train), y = complete_train$y,
        objective = "reg:squarederror", nrounds = 9L, max_depth = 2L,
        learning_rate = 0.35, reg_lambda = 0.7, subsample = 0.8,
        seed = seed, nthreads = 1L
    )
    complete_new <- stats::complete.cases(newdata[, .(x, group)])
    expected <- rep(NA_real_, nrow(newdata))
    expected[complete_new] <- stats::predict(reference, design(newdata[complete_new]))

    set.seed(seed)
    result <- learn_xgboost(
        "y ~ x * group", train, newdata, save_fitted_objects = TRUE,
        xgboost.objective = "reg:squarederror", xgboost.nrounds = 9L,
        xgboost.max_depth = 2L, xgboost.learning_rate = 0.35,
        xgboost.reg_lambda = 0.7, xgboost.subsample = 0.8,
        xgboost.seed = seed, xgboost.nthreads = 1L
    )
    expect_equal(result$predicted_values, expected, tolerance = 1e-6)
    row_order <- c(4, 3, 6, 1, 5, 2)
    reused <- learn_xgboost("y ~ x * group", train, newdata[row_order], reuse_fit = result)
    expect_equal(reused$predicted_values, expected[row_order], tolerance = 1e-6)
})

expect_quiet_learner_error <- function(expr, regexp) {
    printed <- capture.output(expect_error(expr, regexp), type = "message")
    expect_length(printed, 0L)
}

test_that("ranger reports failed fits and incompatible saved model objects", {
    skip_if_not_installed("ranger")
    d <- data.table::data.table(x = 1:8, y = rep(0:1, 4))
    expect_error(learn_ranger("y ~ x", d, d, num.trees = 0L, num.threads = 1L),
                 "Could not fit model with ranger")
    for (type in c("Probability estimation", "Regression")) {
        invalid <- list(fit = list(treetype = type))
        expect_quiet_learner_error(learn_ranger("y ~ x", d, d, reuse_fit = invalid),
                                   "Ranger prediction failed")
    }
})

test_that("xgboost reports errors in fitting, building the design, and prediction", {
    skip_if_not_installed("xgboost", minimum_version = "3.1.2")
    d <- data.table::data.table(x = 1:8, y = rep(0:1, 4))
    expect_error(learn_xgboost("y ~ x", d, d, xgboost.objective = "invalid",
                               xgboost.nthreads = 1L),
                 "Could not fit model with xgboost")
    invalid <- list(fit = list())
    expect_quiet_learner_error(learn_xgboost("y ~ x", d, d, reuse_fit = invalid),
                               "Xgboost prediction failed")
    newdata <- data.table::data.table(x = factor(rep("only level", 4)))
    expect_quiet_learner_error(learn_xgboost("y ~ x", d, newdata, reuse_fit = invalid),
                               "Problem with the design matrix")
})

test_that("tree learners identify the missing optional dependency", {
    ranger <- without_package(learn_ranger, "ranger")
    xgboost <- without_package(learn_xgboost, "xgboost")
    expect_error(ranger(), "Package 'ranger' is required")
    expect_error(xgboost(), "Package 'xgboost' is required")
})
