# Candidate learners returning pre-specified risk scores make the optimal
# weights calculable by hand, independently of a fitting library.
risk_score_learner <- function(character_formula, data, intervened_data, column,
                               save_fitted_objects = FALSE) {
    result <- list(predicted_values = intervened_data[[column]])
    if (save_fitted_objects) result$fit <- column
    result
}

test_that("discrete, NNLS and IPA ensembles return known weights and predictions", {
    d <- data.table::data.table(id = 1:40, a = rep(c(.2, .8), 20),
                                b = rep(c(.4, .6), 20))
    d$y <- (d$a + 2 * d$b) / 3
    newdata <- data.table::data.table(a = c(.15, .85), b = c(.45, .55))
    library <- list(list(name = "a_score", fun = risk_score_learner, args = list(column = "a")),
                    list(name = "b_score", fun = risk_score_learner, args = list(column = "b")))
    # The Brier scores are 4/225 and 1/225. The null score is 1/36.
    weights <- list(discrete = c(a_score = 0, b_score = 1),
                    nnls = c(a_score = 1/3, b_score = 2/3),
                    ipa = c(a_score = .3, b_score = .7))
    for (method in names(weights)) {
        fit <- superlearn(folds = 3, seed = 13, learners = library, parse_learners = FALSE,
                          character_formula = "y ~ a + b", outcome_variable = "y",
                          id_variable = "id", data = d, intervened_data = newdata,
                          ensemble_method = method, save_fitted_objects = TRUE)
        expected <- weights[[method]]
        expect_equal(fit$fit_summary[names(expected)], expected)
        expect_equal(as.numeric(fit$predicted_values),
                     newdata$a * expected[[1]] + newdata$b * expected[[2]])
        expect_named(fit$fit, names(expected)[expected > 0])
    }
})

test_that("IPA falls back to the training mean when both candidates are worse", {
    d <- data.table::data.table(id = 1:40, y = rep(0:1, 20))
    d$wrong <- 1 - d$y
    library <- list(list(name = "first", fun = risk_score_learner, args = list(column = "wrong")),
                    list(name = "second", fun = risk_score_learner, args = list(column = "wrong")))
    fit <- superlearn(folds = 2, seed = 13, learners = library, parse_learners = FALSE,
                      character_formula = "y ~ wrong", outcome_variable = "y",
                      outcome_variable_name = "outcome", id_variable = "id",
                      data = d, intervened_data = d[1:3], ensemble_method = "ipa",
                      diagnostics = list(), save_fitted_objects = TRUE)
    expect_equal(fit$fit_summary, c(first = 0, second = 0))
    expect_equal(fit$predicted_values, rep(0.5, 3))
    expect_equal(fit$diagnostics$superlearner_null_model$Outcome, "outcome")
    expect_null(fit[["fit"]])
})

test_that("super learner records constant fold predictors", {
    d <- data.table::data.table(id = 1:40, y = rep(0:1, 20), x = .5)
    candidates <- list(list(name = "score", fun = risk_score_learner, args = list(column = "x")))
    fit_once <- function(data) {
        superlearn(folds = 2, seed = 13, learners = candidates, parse_learners = FALSE,
                   character_formula = "y ~ x", outcome_variable = "y", outcome_variable_name = "Y",
                   id_variable = "id", data = data, intervened_data = data[1:2], diagnostics = list())
    }
    fit <- fit_once(d)
    expect_equal(fit$diagnostics$constant_predictors$Info, rep("x", 2))
    expect_equal(fit$predicted_values, c(.5, .5))
})

test_that("a constant training outcome is the held-out prediction for every learner", {
    d <- data.table::data.table(id = 1:40, x = 1:40, y = .4)
    set.seed(13)
    fold <- sample(1:2, nrow(d), replace = TRUE)
    d$y[which(fold == 1L)[1:2]] <- c(.2, .6)

    score <- function(intervened_data, value, ...) {
        list(predicted_values = rep(value, nrow(intervened_data)))
    }
    candidates <- list(list(name = "low", fun = score, args = list(value = .2)),
                       list(name = "high", fun = score, args = list(value = .6)))
    fit <- superlearn(folds = 2, seed = 13, learners = candidates, parse_learners = FALSE,
                      character_formula = "y ~ x", outcome_variable = "y", id_variable = "id",
                      data = d, intervened_data = d[1:2], ensemble_method = "nnls")
    # Both held-out folds have mean .4. Training on the constant fold gives
    # predictions (.4, .4); training on the other fold gives (.2, .6).
    # Fitting both fold means therefore requires coefficients (.5, .5).
    expect_equal(fit$fit_summary, c(low = .5, high = .5))
})

test_that("super learner reports character predictors and fold learner failures", {
    d <- data.table::data.table(id = 1:40, y = rep(0:1, 20), x = 1:40)
    failing_learner <- function(...) stop("candidate failed")
    candidates <- list(list(name = "broken", fun = failing_learner))
    fit_once <- function(data) {
        superlearn(folds = 2, seed = 13, learners = candidates, parse_learners = FALSE,
                   character_formula = "y ~ x", outcome_variable = "y", id_variable = "id",
                   data = data, intervened_data = data)
    }
    printed <- capture.output(
        expect_error(fit_once(d), "Learning failed in fold 1 with learner broken: candidate failed"),
        type = "message"
    )
    expect_length(printed, 0L)
    d$x <- as.character(d$x)
    expect_error(fit_once(d), "Offending variables: x")
})

test_that("the discrete learner selects by held-out predictions and refits on all data", {
    d <- data.table::data.table(id = 1:80, x = rep(0:1, each = 40),
                                y = c(rep(c(0, 0, 0, 1), 10), rep(c(0, 1, 1, 1), 10)))
    # This learner is perfect on training IDs, but predicts the training mean
    # for new IDs. Its held-out Brier score is worse than the logistic model's.
    memorize <- function(character_formula, data, intervened_data, ...) {
        idx <- match(intervened_data$id, data$id)
        p <- data$y[idx]
        p[is.na(idx)] <- mean(data$y)
        list(predicted_values = p)
    }
    candidates <- list(list(name = "memorize", fun = memorize),
                       list(name = "glm", fun = learn_glm))
    fit <- superlearn(folds = 4, seed = 13, learners = candidates, parse_learners = FALSE,
                      character_formula = "y ~ x", outcome_variable = "y", id_variable = "id",
                      data = d, intervened_data = data.table::data.table(id = c(81, 82), x = c(1, 0)))
    expect_equal(fit$fit_summary, c(glm = 1, memorize = 0))
    expect_equal(as.numeric(fit$predicted_values), c(3/4, 1/4), tolerance = 1e-6)
})

test_that("factor predictors and named learner specifications work in cross-validation", {
    d <- data.table::data.table(id = 1:120, group = factor(rep(letters[1:3], each = 40)),
                                y = c(rep(c(1, 0, 0, 0), 10), rep(0:1, 20),
                                      rep(c(0, 1, 1, 1), 10)))
    # The two learners use either the three observed group risks or the overall risk.
    candidates <- list(groups = list(fun = "glm"),
                       intercept = list(fun = "glm", learn_variables = "absent"))
    expect_warning(
        fit <- superlearn(folds = 3, seed = 13, learners = candidates,
                          character_formula = "y ~ group", outcome_variable = "y",
                          id_variable = "id", data = d, intervened_data = d[c(120, 60, 1)]),
        "ensemble_method"
    )
    expect_equal(fit$fit_summary, c(groups = 1, intercept = 0))
    expect_equal(as.numeric(fit$predicted_values), c(.75, .5, .25), tolerance = 1e-6)
    expect_error(superlearn(folds = 3, learners = list("glm"), data = d),
                 "length[(]learners[)] > 1")
})
