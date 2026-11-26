test_that("malformed learner specifications identify the invalid field", {
    expect_error(parse_learners(c("glm", "")), "non-empty strings")
    expect_error(parse_learners(list(fun = "glm")), "must be a list that contains 'name' and 'fun'")
    expect_error(parse_learners(list(learners = "glm", folds = 2)), "must be a list of learner")
    expect_error(parse_learners(list(learners = list("glm"), folds = "two")), "'folds' must be numeric")
    expect_error(parse_learners(42), "Invalid 'learners' specification")
    expect_error(parse_learners(list(name = "bad", fun = 1)), "must specify a learner function")
    for (fun in list(NA_character_, "", c("glm", "glmnet"))) {
        expect_error(parse_learners(list(name = "bad", fun = fun)), "single non-empty")
    }
})

test_that("learner names, arguments and cross-validation defaults are preserved", {
    single <- parse_learners(list(name = NULL, fun = learn_glm, maxit = 12))
    expect_identical(single, list(name = "unnamed", fun = learn_glm, args = list(maxit = 12)))
    candidates <- list(list(fun = "glm"), list(name = "forest", fun = "ranger", num.trees = 5L))
    expect_warning(
        ensemble <- parse_learners(list(learners = candidates, folds = numeric(), ensemble_method = "ipa")),
        "Defaults to 10-fold"
    )
    expect_identical(ensemble$learners[[1]]$name, "noname_learner_1")
    expect_identical(ensemble$learners[[2]]$args, list(num.trees = 5L))
    expect_equal(ensemble$args$folds, 10)
    expect_identical(ensemble$args$ensemble_method, "ipa")
    expect_warning(
        named <- parse_learners(list(learners = list(logistic = list(fun = "glm")), folds = 2L)),
        "Defaults to .discrete."
    )
    expect_identical(named$learners$logistic$name, "logistic")
    expect_identical(named$args$ensemble_method, "discrete")
})

test_that("character vectors select a cross-validated library of learners", {
    parsed <- parse_learners(c("glm", "ranger"))
    expect_identical(parsed$name, "superlearn")
    expect_equal(parsed$args$folds, 10)
    expect_identical(vapply(parsed$learners, "[[", "", "fun"), c("learn_glm", "learn_ranger"))
})
