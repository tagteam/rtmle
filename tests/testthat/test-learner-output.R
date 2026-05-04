library(testthat)
library(rtmle)
library(data.table)

test_that("learners return list outputs without prediction attributes",{
    d <- data.table(Y = rep(c(0, 1), 10),
                    A = rep(c(0, 1, 1, 0), 5),
                    L = seq(-1, 1, length.out = 20))

    glm_out <- learn_glm("Y ~ A + L", data = d, intervened_data = d)
    expect_type(glm_out,"list")
    expect_identical(names(glm_out)[1],"predicted_values")
    expect_null(attributes(glm_out$predicted_values))
    expect_equal(length(glm_out$predicted_values),NROW(d))
    expect_false(is.null(glm_out$fit))
    expect_false(is.null(glm_out$object))

    glmnet_out <- learn_glmnet("Y ~ A + L", data = d, intervened_data = d,
                               lambda = 0.01)
    expect_type(glmnet_out,"list")
    expect_identical(names(glmnet_out)[1],"predicted_values")
    expect_null(attributes(glmnet_out$predicted_values))
    expect_equal(length(glmnet_out$predicted_values),NROW(d))
    expect_false(is.null(glmnet_out$fit))
    expect_false(is.null(glmnet_out$object))
    expect_false(is.null(glmnet_out$selected.lambda))
})

test_that("run_rtmle can store fitted learner objects",{
    set.seed(8)
    n <- 80
    W <- rnorm(n)
    A <- rbinom(n,1,plogis(-1 + 2 * W))
    Y <- rbinom(n,1,plogis(W + A))

    make_x <- function(){
        x <- rtmle_init(time_grid = 0:1,
                        name_id = "id",
                        name_outcome = "Y",
                        name_competing = NULL,
                        name_censoring = NULL)
        x$prepared_data <- data.table(id = 1:n,W = W,A_0 = A,Y_1 = Y)
        x$names$name_baseline_covariates <- "W"
        x$names$name_time_covariates <- "A"
        x <- protocol(x,name = "A",treatment_variables = "A",intervention = 1)
        x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "A")
        model_formula(x)
    }

    x_summary <- run_rtmle(make_x(),learner = "learn_glm",time_horizon = 1,
                           refit = TRUE,verbose = FALSE)
    summary_fit <- x_summary$models$time_0$A$A_0$fit
    expect_true(is.matrix(summary_fit) || is.data.frame(summary_fit))

    x_object <- run_rtmle(make_x(),learner = "learn_glm",time_horizon = 1,
                          refit = TRUE,verbose = FALSE,
                          save_fitted_objects = TRUE)
    object_fit <- x_object$models$time_0$A$A_0$fit
    expect_true(inherits(object_fit,"glm") || inherits(object_fit,"speedglm"))
})
