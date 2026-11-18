test_that("super learning retains the selected learner fit when requested", {
    saved_flags <- logical()
    learner <- function(character_formula, data, intervened_data, shift = 0,
                         save_fitted_objects = FALSE, ...) {
        saved_flags <<- c(saved_flags, save_fitted_objects)
        fit <- structure(list(shift = shift), class = "test_ensemble_member")
        output <- list(predicted_values = stats::plogis(4 * intervened_data$x + shift),
                       fit_summary = shift)
        if (save_fitted_objects) output$fit <- fit
        output
    }
    candidates <- parse_learners(list(folds = 2, ensemble_method = "discrete",
                                      learners = list(good = list(fun = learner),
                                                      bad = list(fun = learner, shift = 5))))
    d <- data.table::data.table(Y = rep(c(0, 1), 30), x = rep(c(-1, 1), 30))
    run <- function(save) superlearn(
        folds = 2, seed = 15, learners = candidates$learners,
        parse_learners = FALSE, character_formula = "Y ~ x",
        outcome_variable = "Y", id_variable = "id", data = d,
        intervened_data = d, save_fitted_objects = save)
    saved <- run(TRUE)
    expect_named(saved$fit, "good")
    expect_s3_class(saved$fit$good, "test_ensemble_member")
    expect_true(tail(saved_flags, 1L))
    unsaved <- run(FALSE)
    expect_false("fit" %in% names(unsaved))
    expect_equal(unsaved$predicted_values, saved$predicted_values)
    expect_identical(unsaved$fit_summary, saved$fit_summary)
})

test_that("retained glmnet ensemble members are actual fitted models", {
    d <- data.table::data.table(Y = rep(c(0, 1, 0, 1), 20),
                                x = rep(c(-1, 0.5, -0.5, 1), 20),
                                z = seq(-1, 1, length.out = 80))
    candidates <- parse_learners(list(folds = 2, ensemble_method = "discrete",
                                      learners = list(
                                          small = list(fun = "glmnet", lambda = 0.01),
                                          large = list(fun = "glmnet", lambda = 0.1))))
    result <- superlearn(folds = 2, seed = 9, learners = candidates$learners,
                          parse_learners = FALSE, character_formula = "Y ~ x + z",
                          outcome_variable = "Y", id_variable = "id", data = d,
                          intervened_data = d, save_fitted_objects = TRUE)
    active <- names(result$fit_summary)[result$fit_summary > 0]
    expect_named(result$fit, active)
    expect_true(all(vapply(result$fit, inherits, logical(1), "glmnet")))
    expect_true(all(is.finite(result$predicted_values)))
})
