test_that("dynamic rules distinguish fixed actions, natural treatment and unsupported flips", {
    x <- rtmle_init(c(0, 2), name_id = "id", name_time = "time",
                    name_outcome = "Y", continuous_time = TRUE)
    rule <- function(data, intervention_table, time_node) {
        data$A_0[data$L_0 == 1] <- 1
        data
    }
    x <- regime(x, "rule", intervention = 1, intervene_function = rule)
    history <- data.table::data.table(id = 1:4, A_0 = c(0, 0, 1, 1), L_0 = c(0, 1, 0, 1))
    policy <- rtmle:::continuous_intervention_policy(history, x$regimes$rule, 0)
    expect_equal(policy$value, c(0, 1, 1, 1))
    expect_equal(policy$natural, c(TRUE, FALSE, TRUE, FALSE))

    flip <- function(data, intervention_table, time_node) {
        data$A_0 <- 1 - data$A_0
        data
    }
    x <- regime(x, "flip", intervention = 1, intervene_function = flip)
    expect_error(rtmle:::continuous_intervention_policy(history, x$regimes$flip, 0),
                 "modified-treatment-policy")
})

test_that("natural treatment has probability one and is excluded from adherence fitting", {
    d <- continuous_example()
    d$baseline_data$L_0 <- c(0, 1, 0, 1)
    prepared <- prepare_continuous_data(d, 2, c("A", "L"), c("A_0", "L_0"))
    x <- rtmle_init(c(0, 2), name_id = "id", name_time = "time",
                    name_outcome = "Y", continuous_time = TRUE)
    rule <- function(data, intervention_table, time_node) {
        data$A_0[data$L_0 == 1] <- 1
        data
    }
    x <- regime(x, "rule", intervention = 1, intervene_function = rule)
    nuisance <- rtmle:::continuous_propensity(prepared, learner = "learn_glm",
                                              regime = x$regimes$rule)
    # One of the two constrained subjects adheres. Including the two natural
    # subjects in the same intercept model would instead produce 3/4.
    expect_equal(nuisance$data$propensity_0, c(1, 0.5, 1, 0.5))
    x <- regime(x, "rule", intervention = 1, intervene_function = rule,
                adherence_model_strata = "L")
    nuisance <- rtmle:::continuous_propensity(prepared, learner = "learn_glm",
                                              regime = x$regimes$rule)
    expect_equal(nuisance$data$propensity_0, c(1, 0.5, 1, 0.5))
})

test_that("unsuffixed adherence strata use the current pre-decision covariate", {
    baseline <- data.frame(id = 1:8, A_0 = 0, L_0 = rep(0:1, each = 4))
    assignment <- c(1, 0, 0, 0, 1, 1, 1, 0)
    measured <- rep(0:1, 4)
    events <- rbind(data.frame(id = 1:8, time = 1, event = "A", A = assignment, L = measured),
                    data.frame(id = 1:8, time = 2, event = c("Y", rep("tauend", 7)),
                               A = assignment, L = measured))
    prepared <- prepare_continuous_data(list(baseline_data = baseline, timevarying_data = events),
                                        2, c("A", "L"), c("A_0", "L_0"), min_events = 0)
    x <- rtmle_init(c(0, 2), name_id = "id", name_time = "time",
                    name_outcome = "Y", continuous_time = TRUE)
    x <- regime(x, "rule", 1, intervene_function = intervene, adherence_model_strata = "L")
    mean_learner <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = rep(mean(data$rtmle_outcome), nrow(intervened_data)))
    }
    fit <- rtmle:::continuous_propensity(prepared, regime = x$regimes$rule, learner = mean_learner)
    expect_equal(fit$data$propensity_1, rep(c(0.75, 0.25), 4))
    expect_equal(rtmle:::continuous_adherence_strata("L_1", 1L, names(prepared$wide_data),
        prepared$baseline_covariates, prepared$time_covariates), "L_1")
})

test_that("continuous adherence strata fit separate propensities", {
    d <- continuous_example()
    d$baseline_data <- data.frame(id = 1:8, A_0 = c(1, 0, 0, 0, 1, 1, 1, 0),
                                  L_0 = rep(0:1, each = 4))
    d$timevarying_data <- data.frame(id = 1:8, time = 2, event = "tauend",
                                     A = d$baseline_data$A_0, L = d$baseline_data$L_0)
    prepared <- prepare_continuous_data(d, 2, c("A", "L"), c("A_0", "L_0"))
    x <- rtmle_init(c(0, 2), name_id = "id", name_time = "time",
                    name_outcome = "Y", continuous_time = TRUE)
    x <- regime(x, "rule", intervention = 1, intervene_function = intervene,
                adherence_model_strata = "L")
    fit <- rtmle:::continuous_propensity(prepared, regime = x$regimes$rule)
    expect_equal(fit$data$propensity_0, rep(c(0.25, 0.75), each = 4))
    expect_length(fit$models$treatment$A_0$fit_summary, 2)
})

test_that("continuous adherence strata require observed pre-decision variables", {
    x <- rtmle_init(c(0, 2), name_id = "id", name_time = "time",
                    name_outcome = "Y", continuous_time = TRUE)
    expect_error(regime(x, "bad", 1, adherence_model_strata = "L"), "only for dynamic")
    expect_error(regime(x, "bad", 1, intervene_function = intervene,
                        adherence_model_strata = NA_character_), "character vector")
    prepared <- prepare_continuous_data(continuous_example(), 2,
                                         c("A", "L"), c("A_0", "L_0"))
    make_rule <- function(strata) regime(x, "rule", 1, intervene_function = intervene,
                                         adherence_model_strata = strata)$regimes$rule
    expect_error(rtmle:::continuous_propensity(prepared, regime = make_rule("unknown")), "Unknown")
    expect_error(rtmle:::continuous_propensity(prepared, regime = make_rule("L_1")), "pre-decision")
    expect_error(rtmle:::continuous_propensity(prepared, regime = make_rule("A")), "baseline treatment")
    prepared$wide_data$L_0[1] <- NA
    expect_error(rtmle:::continuous_propensity(prepared, regime = make_rule("L")), "must be observed")
})


test_that("dynamic rules retain numerically suffixed fixed baseline predictors", {
    d <- continuous_example()
    d$baseline_data$year_2024 <- c(0, 1, 0, 1)
    prepared <- prepare_continuous_data(d, 2, c("A", "L"), c("A_0", "L_0", "year_2024"))
    x <- rtmle_init(c(0, 2), name_id = "id", name_time = "time", name_outcome = "Y",
                    name_competing = "D", name_censoring = "C", continuous_time = TRUE)
    rule <- function(data, intervention_table, time_node) {
        data$A_0[data$year_2024 == 1] <- 1
        data
    }
    x <- regime(x, "rule", intervention = 1, intervene_function = rule)
    nuisance <- rtmle:::continuous_propensity(prepared, learner = "learn_glm", regime = x$regimes$rule)
    expect_equal(nuisance$data$propensity_0, c(1, 0.5, 1, 0.5))

    x <- continuous_add_example_data(x, d)
    x <- prepare_rtmle_data(x)
    x <- target(x, "risk", regimes = "rule")
    x <- model_formula(x, verbose = FALSE)
    expect_error(run_rtmle(x, estimator = "one-step", intervene_all_sequential_regression = TRUE),
                 "not supported for dynamic")
    result <- run_rtmle(x, estimator = "one-step")
    expect_true(is.finite(result$estimate$Main_analysis$Estimate))
})
