library(testthat)

test_that("learner method names are expanded to learn_ functions", {
    expect_identical(
        parse_learners("glm")[c("name", "fun")],
        list(name = "learn_glm", fun = "learn_glm")
    )
    expect_identical(
        parse_learners("learn_glm")[c("name", "fun")],
        list(name = "learn_glm", fun = "learn_glm")
    )

    parsed <- parse_learners(list(name = "mylearner", fun = "glm"))
    expect_identical(parsed$name, "mylearner")
    expect_identical(parsed$fun, "learn_glm")

    for (method in c("glm", "glmnet", "ranger", "xgboost")) {
        expect_identical(
            parse_learners(method)$fun,
            paste0("learn_", method)
        )
    }
})


test_that("missing learner functions report the normalized function name", {
    expected_message <- paste0(
        "cannot find function learn_yyy",
        " maybe you have to write this first?"
    )

    expect_error(parse_learners("yyy"), expected_message, fixed = TRUE)
    expect_error(parse_learners("learn_yyy"), expected_message, fixed = TRUE)
    expect_error(
        parse_learners(list(name = "mylearner", fun = "yyy")),
        expected_message,
        fixed = TRUE
    )
})
