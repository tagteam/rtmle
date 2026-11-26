test_that("fitter records why it uses an empirical probability", {
    d <- data.table::data.table(Y = c(0, 0, 0, 1), A = 1L)
    fit <- function(formula, threshold, diagnostics = NULL) {
        rtmle:::fitter(1, parse_learners("glm"), formula, d, d[1:2],
                       "id", threshold, seed = 1, diagnostics = diagnostics)
    }
    first <- fit("Y ~ A", 2)
    expect_equal(first$predicted_values, c(.25, .25))
    expect_equal(first$diagnostics$constant_outcome_variables$Minority_count, 1)
    second <- fit("Y ~ A", 2, first$diagnostics)
    expect_equal(nrow(second$diagnostics$constant_outcome_variables), 2)

    # Here the outcome varies enough to fit, but no predictor varies.
    empty <- fit("Y ~ 1", 1)
    expect_equal(empty$fit_summary, .25)
    expect_equal(empty$diagnostics$empty_rhs$Event, "Empty rhs of formula")
    with_constant <- fit("Y ~ A", 1, empty$diagnostics)
    expect_equal(with_constant$diagnostics$constant_predictor_variables$Info, "A")
    expect_equal(nrow(with_constant$diagnostics$empty_rhs), 2)
    expect_equal(fit("Y ~ A", 1)$predicted_values, c(.25, .25))

    d$Y <- 0
    expect_equal(fit("Y ~ A", 1)$diagnostics$constant_outcome_variables$Minority_count, 0)
})

test_that("fitter returns learner warnings and identifies failed outcome models", {
    d <- data.table::data.table(Y = rep(0:1, 10), x = 1:20)
    d$rtmle_predicted_outcome <- d$Y
    apply_learner <- function(learner, outcome = "Y", data = d) {
        rtmle:::fitter(3L, learner, paste(outcome, "~ x"), data, data,
                       "id", minority_threshold = 1, seed = 8)
    }
    noisy <- function(data, intervened_data, ...) {
        warning("few observations")
        warning("fit is singular")
        list(predicted_values = rep(mean(data$rtmle_outcome), nrow(intervened_data)))
    }
    output <- apply_learner(list(name = "custom", fun = noisy))
    expect_equal(output$predicted_values, rep(.5, 20))
    expect_match(output$diagnostics$learner_warnings$warning,
                 "few observations.*fit is singular")
    failing <- list(name = "custom", fun = function(...) stop("solver failure"))
    ensemble <- parse_learners(list(folds = 2, ensemble_method = "discrete",
                                    learners = list("glm")))
    bad <- data.table::copy(d)
    bad$x <- as.character(bad$x)
    for (outcome in c("Y", "rtmle_predicted_outcome")) {
        label <- if (outcome == "Y") "outcome is: Y" else "intervention node \\(time\\): 3"
        expect_error(apply_learner(failing, outcome),
                     paste0("Failed to apply.*\n.*", label, "\nsolver failure"))
        expect_error(apply_learner(ensemble, outcome, bad),
                     paste0("Failed to superlearn/crossfit.*\n.*", label, "\n.*\nOffending variables: x"))
    }
})

test_that("formula helpers replace responses and remove terms with recorded dimensions", {
    renamed <- rtmle:::change_outcome_in_formula(Y ~ age + I(age^2), "Z")
    expect_equal(attr(renamed, "all_vars"), c("Z", "age"))
    expect_equal(as.character(renamed), "Z ~ age + I(age^2)")
    reduced <- rtmle:::delete_variables_from_formula(Y ~ age + x, "x")
    expect_equal(as.character(reduced), "Y ~ age")
    expect_equal(attr(reduced, "number_rhs_variables"), 1L)
    empty <- rtmle:::delete_variables_from_formula("Y ~ x", "x")
    expect_equal(as.character(empty), "Y ~ 1")
    expect_equal(attr(empty, "number_rhs_variables"), 0L)
    expect_error(rtmle:::change_outcome_in_formula(c("Y ~ A", "Z ~ B"), "C"),
                 "single character string")
    expect_error(rtmle:::delete_variables_from_formula(c("Y ~ A", "Z ~ B"), "A"),
                 "single character string")
    expect_equal(make_exclusion_rule(1, "Y", c("A", "B"), recursive = TRUE),
                 list(Y_0 = c("A_0", "B_0"), Y_1 = c("A_0", "A_1", "B_0", "B_1")))
    expect_equal(make_exclusion_rule(1, "Y", "A"), list(Y_0 = "A_0", Y_1 = "A_1"))
})

test_that("history restrictions specify the intended propensity regressions", {
    args <- list(timepoint = 2, available_names = c("L_0", "L_1", "L_2", "age"),
                 name_outcome_variable = "A_2", outcome_value = 1,
                 name_baseline_covariates = "age", name_time_covariates = "L",
                 Markov = "L", handle_concomitant_variables = "joint")
    form <- function(...) do.call(rtmle:::formalize, c(args, list(...)))$A_2$formula
    expect_equal(form(), "I(1*(A_2==1)) ~ age + L_1")
    # If current L is known at assignment, the Markov model replaces lagged L.
    expect_equal(form(additional_variables = "L_2"), "I(1*(A_2==1)) ~ age + L_2")
    expect_equal(form(exclusion_rules = list(A_2 = "age"),
                      inclusion_rules = list(A_2 = "L_0")), "I(1*(A_2==1)) ~ L_1 + L_0")
    expect_warning(form(inclusion_rules = list(A_2 = "missing")), "did not match any variable")

    joint <- rtmle:::formalize(0, character(), c("A_0", "B_0"), c(1, 0),
                               character(), character(),
                               handle_concomitant_variables = "sequential")
    expect_equal(joint[[1]]$A_0$formula, "I(1*(A_0==1)) ~ 1")
    expect_equal(joint[[2]]$B_0$formula, "I(1*(B_0==0)) ~ A_0")
})

test_that("model_formula diagnoses invalid regimes, histories and formula names", {
    x <- rtmle_init(0:2, name_id = "id", name_outcome = "Y")
    expect_error(model_formula(x), "No regimes registered")
    data(rtmle_object)
    expect_error(model_formula(rtmle_object, Markov = "unmeasured"),
                 "do not match time_covariates")
    # An imported baseline column name containing a formula separator produces
    # an ambiguous response when a dynamic propensity formula is constructed.
    baseline <- data.table::copy(rtmle_object$data$baseline_data)
    baseline[["baseline ~ age"]] <- baseline$age
    x <- add_baseline_data(rtmle_object, baseline)
    expect_warning(x <- prepare_rtmle_data(x, verbose = FALSE), "Overwriting existing")
    for (method in c("independent", "sequential")) {
        y <- regime(x, name = "AB", treatment_variables = c("A", "B"),
                    intervention = c(1, 0), multiple_treatment_factorization = method,
                    intervene_function = treat_unless(contra_indication = "bleeding",
                                                     action = c(0, NA)),
                    verbose = FALSE)
        expect_error(model_formula(y, verbose = FALSE),
                     "Could not identify the response in propensity formula")
    }
})

test_that("simulation models use specified distributions and require event scales", {
    variables <- list(binary = "binary", x = "normal", positive = "lognormal",
                      fraction = "uniform", wait = "exponential", time = "weibull",
                      fixed = "constant")
    parameters <- list(intercept_binary = 1, intercept_positive = .2,
                       var_positive = .5, scale_wait = .4, scale_time = .3,
                       intercept_fixed = 3, effect_x_positive = 1)
    model <- rtmle:::make_regression_model(variables, parameters)
    d <- lava::sim(model, n = 2500, seed = 1)
    expect_true(all(d$binary %in% 0:1))
    expect_equal(unique(d$fixed), 3)
    expect_true(all(d$fraction >= 0 & d$fraction <= 1))
    expect_true(all(d$positive > 0 & d$wait > 0 & d$time > 0))
    residual <- log(d$positive) - d$x
    expect_lt(abs(mean(residual) - .2), .05)
    expect_lt(abs(stats::sd(residual) - .5), .05)
    expect_error(rtmle:::make_regression_model(list(time = "weibull"), list()),
                 "missing parameter scale_time")
})
