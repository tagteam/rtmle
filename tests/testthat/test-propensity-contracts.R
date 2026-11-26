make_propensity_fixture <- function(factorization, rule = NULL, strata = NULL,
                                    baseline_data = NULL) {
    x <- rtmle_init(time_grid = 0:1, name_id = "id", name_outcome = "Y",
                    minority_threshold = 1)
    a <- rep(0:1, each = 8)
    b <- c(rep(0:1, c(6, 2)), rep(0:1, c(2, 6)))
    x <- add_wide_data(x, outcome_data = data.frame(id = 1:16, Y_1 = rep(0:1, 8)),
                       timevar_data = list(A = data.frame(id = 1:16, A_0 = factor(a)),
                                           B = data.frame(id = 1:16, B_0 = factor(b))))
    if (!is.null(baseline_data)) x <- add_baseline_data(x, baseline_data)
    x <- prepare_rtmle_data(x, verbose = FALSE)
    x <- regime(x, "both", intervention = c(1,1), treatment_variables = c("A","B"),
                multiple_treatment_factorization = factorization,
                intervene_function = rule, adherence_model_strata = strata)
    model_formula(x, verbose = FALSE)
}
fit_baseline_propensity <- function(x, learner = parse_learners("learn_glm"), progressbar = 0) {
    rtmle:::intervention_probabilities(x, "both", max_intervention_node = 0L,
                               refit = TRUE, learner = learner, seed = 17L,
                               progressbar = progressbar, save_fitted_objects = TRUE)
}

test_that("joint, sequential and independent propensity models give known probabilities", {
    # Of 16 participants, eight take A, eight take B, and six take both.
    # Thus P(A,B) = 6/16, P(B|A) = 6/8, and P(A)P(B) = 1/4.
    expected <- list(joint = 6/16, sequential = c(1/2, 6/16),
                     independent = c(1/2, 1/4))
    for (factorization in names(expected)) {
        x <- make_propensity_fixture(factorization)
        fitted <- fit_baseline_propensity(x)
        probs <- fitted$regimes$both$cumulative_intervention_probs
        oracle <- matrix(rep(expected[[factorization]], each = 16), nrow = 16)
        expect_equal(unname(probs), oracle, tolerance = 1e-7)
        expect_equal(sum(fitted$regimes$both$intervention_match[, 1]), 6)
    }
})

test_that("dynamic adherence probabilities are specific to baseline strata", {
    static_rule <- function(data, intervention_table, time_node) {
        intervene(data, intervention_table, time_node)
    }
    baseline <- data.frame(id = 1:16, W = rep(c(0, 0, 1, 1), 4))
    x <- make_propensity_fixture("joint", static_rule, "W", baseline)
    x$followup <- NULL
    fitted <- fit_baseline_propensity(x)
    # Two participants with W = 0 and four with W = 1 take both treatments.
    expect_equal(as.numeric(fitted$regimes$both$intervention_probs),
                 ifelse(x$prepared_data$W == 0, 2/8, 4/8))
    fit_summary <- fitted$models$time_0$both[["A_0,B_0"]]$fit_summary
    expect_s3_class(fit_summary, "stratified_adherence_probability")
    expect_named(fit_summary, c("0", "1"))
    expect_equal(unlist(unclass(fit_summary), use.names = FALSE), c(2/8, 4/8))

    x$prepared_data$W[1] <- NA
    expect_error(fit_baseline_propensity(x), "strata.*must be observed")
})

test_that("propensity fitting diagnoses invalid model specifications and predictions", {
    x <- make_propensity_fixture("joint")
    x$regimes$both$intervene_function <- NULL
    expect_error(fit_baseline_propensity(x), "No intervene function defined for regime both")

    # A joint-treatment model uses two treatment columns. Both must be present.
    x <- make_propensity_fixture("joint", intervene)
    names(x$models$time_0$both) <- "A_0,missing"
    expect_error(fit_baseline_propensity(x), "Cannot construct adherence.*missing")
    names(x$models$time_0$both) <- "A_0,B_0"
    x$models$time_0$both[[1]]$formula <- "A_0"
    expect_error(fit_baseline_propensity(x), "identify the response in propensity formula")

    # Check the prediction contract using a learner with prescribed outputs.
    x <- make_propensity_fixture("sequential")
    fixed <- function(character_formula, data, intervened_data, probabilities, ...) {
        list(predicted_values = rep_len(probabilities, nrow(intervened_data)))
    }
    fitted_with <- function(probabilities) {
        fit_baseline_propensity(x, list(name = "fixed", fun = fixed,
                                       args = list(probabilities = probabilities)))
    }
    expect_error(fitted_with(NA_real_), "returned missing values")
    expect_error(fitted_with(0), "predictions are exactly zero")
    bounded <- fitted_with(c(-0.1, 1.1))
    expect_equal(as.numeric(bounded$regimes$both$intervention_probs[, "B_0"]),
                 rep(0:1, 8))
    expect_equal(bounded$diagnostics$probabilities_off_range$variable, "B_0")
})

test_that("a treatment task with no one at risk has no fitted model", {
    x <- make_propensity_fixture("joint")
    x$followup$last_interval[] <- -1L
    fitted <- fit_baseline_propensity(x)
    expect_s3_class(fitted$models$time_0$both[[1]]$fit_summary, "no_propensity_model")
    expect_true(all(is.na(fitted$regimes$both$intervention_probs)))
})

test_that("propensity checks detect incompatible levels and duplicate sequential steps", {
    x <- make_propensity_fixture("joint")
    # The desired level 1 cannot be assigned if the data retain only level 0.
    x$prepared_data$A_0 <- factor(rep(0, nrow(x$prepared_data)))
    expect_error(fit_baseline_propensity(x), "changes a nominal treatment value")
    x <- make_propensity_fixture("sequential")
    steps <- x$models$time_0$both
    x$models$time_0$both <- c(steps, steps[1])
    expect_error(fit_baseline_propensity(x), "locate the sequential propensity model")
})

test_that("the progress display identifies the regime being fitted", {
    x <- make_propensity_fixture("joint")
    printed <- capture.output(
        expect_message(fit_baseline_propensity(x, progressbar = 3),
                       "Fitting propensity score and censoring models: both")
    )
    expect_true(any(grepl("100%", printed, fixed = TRUE)))
})
