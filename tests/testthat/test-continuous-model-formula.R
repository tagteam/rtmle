# Every person has the same event ordering, but treatment and covariate
# histories vary. Some people are censored at the last event.
continuous_formula_fixture <- function() {
    n <- 24L
    id <- seq_len(n)
    baseline_A <- rep(c(0, 0, 1, 1, 0, 1), length.out = n)
    baseline_L <- id %% 4
    first_A <- rep(c(0, 1, 1, 0, 0, 1), length.out = n)
    last_A <- rep(c(0, 1, 1, 0), length.out = n)
    L1 <- (2 * id + id %% 3) %% 7
    L2 <- (3 * id + id %% 4) %% 11
    L3 <- (id %% 5 + id %% 3) %% 7
    terminal <- rep(c("Y", "C", "tauend"), length.out = n)

    x <- rtmle_init(c(0, 4), name_id = "id", name_time = "time",
                    name_outcome = "Y", name_competing = "D",
                    name_censoring = "C", continuous_time = TRUE)
    x <- add_baseline_data(x, data.frame(id = id, A = baseline_A,
                                          L = baseline_L, age = 40 + id,
                                          offset_value = 1 + id / 10,
                                          group = factor(rep(c("a", "a", "b", "b"), length.out = n))))
    events <- data.frame(
        id = rep(id, each = 4), time = rep(1:4, n),
        event = as.vector(rbind(rep("A", n), rep("visit", n), rep("A", n), terminal)),
        A = as.vector(rbind(first_A, first_A, last_A, last_A)),
        L = as.vector(rbind(L1, L2, L3, L3))
    )
    x <- add_continuous_data(x, event_data = events, treatment_events = "A",
                             covariate_events = "visit", time_covariates = "L")
    x <- prepare_rtmle_data(x, min_events = 0)
    x <- regime(x, "Always_A", intervention = 1, treatment_variables = "A")
    target(x, "risk", regimes = "Always_A")
}

# Capture formulas actually passed to nuisance learners. Returning constant,
# valid predictions keeps the test focused on formula construction.
continuous_formula_spy <- function() {
    captured <- list()
    learner <- function(character_formula, data, intervened_data,
                        time_variable = NULL, ...) {
        captured[[length(captured) + 1L]] <<- list(
            formula = character_formula, time_variable = time_variable)
        list(predicted_values = rep(0.5, nrow(intervened_data)), fit = list())
    }
    list(learner = learner, calls = function() captured)
}

continuous_formula_fit <- function(x, treatment, censor, refit = TRUE) {
    run_rtmle(x, learner = "learn_glm", estimator = "one-step",
              propensity_learner = treatment$learner,
              censoring_learner = list(name = "censor_spy", fun = censor$learner,
                                       marginal_censoring = FALSE),
              refit = refit)
}

continuous_formula_fit_marginal <- function(x, treatment, censor, formula = NULL) {
    censoring <- list(name = "censor_spy", fun = censor$learner,
                      marginal_censoring = TRUE)
    if (!is.null(formula)) censoring$marginal_censoring_formula <- formula
    run_rtmle(x, learner = "learn_glm", estimator = "one-step",
              propensity_learner = treatment$learner, censoring_learner = censoring)
}

continuous_formula_rhs <- function(formula) {
    all.vars(stats::as.formula(formula)[[3L]])
}

continuous_formula_outcome <- function(fit, event) {
    fit$sequential_outcome_regression$risk$Always_A$time_horizon_1[[event]]$formula
}

continuous_formula_censor <- function(spy, event) {
    calls <- spy$calls()
    wanted <- Filter(function(call) identical(call$time_variable, paste0("time_", event)), calls)
    if (!length(wanted)) stop("No captured censoring model for event ", event)
    wanted[[length(wanted)]]$formula
}

test_that("Markov limits each model to the most recent covariate value", {
    x <- continuous_formula_fixture()
    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    default <- continuous_formula_fit(model_formula(x, verbose = FALSE), treatment, censor)
    outcome_default <- continuous_formula_rhs(continuous_formula_outcome(default, "event_4"))
    expect_true(all(paste0("L_", 0:3) %in% outcome_default))

    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    markov <- continuous_formula_fit(model_formula(x, Markov = "L", verbose = FALSE),
                                      treatment, censor)
    treatment_calls <- treatment$calls()
    expect_gte(length(treatment_calls), 3L)
    baseline_rhs <- continuous_formula_rhs(treatment_calls[[1L]]$formula)
    treatment_rhs <- continuous_formula_rhs(tail(treatment_calls, 1L)[[1L]]$formula)
    outcome_rhs <- continuous_formula_rhs(continuous_formula_outcome(markov, "event_4"))
    censor_rhs <- continuous_formula_rhs(continuous_formula_censor(censor, 4L))
    expect_true("L_0" %in% baseline_rhs)
    for (rhs in list(treatment_rhs, outcome_rhs, censor_rhs)) {
        expect_true("L_3" %in% rhs)
        expect_false(any(paste0("L_", 0:2) %in% rhs))
    }
})

test_that("named exclusion and inclusion rules reach treatment, outcome and censoring models", {
    x <- continuous_formula_fixture()
    rules <- list("^A_3$" = "^L_3$", "^Y_4$" = "^L_2$", "^C_4$" = "^L_2$")
    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    fit <- continuous_formula_fit(model_formula(x, exclusion_rules = rules, verbose = FALSE),
                                  treatment, censor)
    treatment_rhs <- continuous_formula_rhs(tail(treatment$calls(), 1L)[[1L]]$formula)
    outcome_rhs <- continuous_formula_rhs(continuous_formula_outcome(fit, "event_4"))
    censor_rhs <- continuous_formula_rhs(continuous_formula_censor(censor, 4L))
    expect_false("L_3" %in% treatment_rhs)
    expect_false("L_2" %in% outcome_rhs)
    expect_false("L_2" %in% censor_rhs)
    expect_true("L_1" %in% treatment_rhs)
    expect_true("L_1" %in% outcome_rhs)
    expect_true("L_1" %in% censor_rhs)

    # As in the discrete formula workflow, inclusion is applied after exclusion.
    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    x <- model_formula(x, exclusion_rules = rules,
                       inclusion_rules = list("^Y_4$" = "L_2"), verbose = FALSE)
    fit <- continuous_formula_fit(x, treatment, censor)
    expect_true("L_2" %in% continuous_formula_rhs(continuous_formula_outcome(fit, "event_4")))

    # An inclusion rule can deliberately add an older value after Markov
    # has selected the most recent one.
    x <- model_formula(x, Markov = "L", inclusion_rules = list("^Y_4$" = "^L_1$"), verbose = FALSE)
    fit <- continuous_formula_fit(x, continuous_formula_spy(), continuous_formula_spy())
    rhs <- continuous_formula_rhs(continuous_formula_outcome(fit, "event_4"))
    expect_true(all(c("L_1", "L_3") %in% rhs))
    expect_false(any(c("L_0", "L_2") %in% rhs))
})

test_that("changing model formulas invalidates cached fits even with refit FALSE", {
    x <- model_formula(continuous_formula_fixture(), verbose = FALSE)
    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    first <- continuous_formula_fit(x, treatment, censor, refit = FALSE)
    first_n <- length(treatment$calls())
    expect_true("L_2" %in% continuous_formula_rhs(continuous_formula_outcome(first, "event_4")))

    changed <- model_formula(first, Markov = "L", verbose = FALSE)
    second <- continuous_formula_fit(changed, treatment, censor, refit = FALSE)
    expect_gt(length(treatment$calls()), first_n)
    expect_false("L_2" %in% continuous_formula_rhs(continuous_formula_outcome(second, "event_4")))
    expect_true("L_3" %in% continuous_formula_rhs(continuous_formula_outcome(second, "event_4")))
})

test_that("formula options reject invalid names and future-only inclusion", {
    x <- continuous_formula_fixture()
    expect_error(model_formula(x, Markov = "unknown", verbose = FALSE), "Markov")
    expect_error(model_formula(x, exclusion_rules = list("L_3"), verbose = FALSE), "[Nn]am|rule")
    expect_error(model_formula(x, inclusion_rules = list("L_3"), verbose = FALSE), "[Nn]am|rule")
    expect_error({
        bad <- model_formula(x, inclusion_rules = list("^A_1$" = "L_3"), verbose = FALSE)
        continuous_formula_fit(bad, continuous_formula_spy(), continuous_formula_spy())
    }, "[Ff]uture|[Uu]navailable|[Hh]istory")
})

test_that("global exclusions cover baseline and indexed time-varying covariates", {
    x <- model_formula(continuous_formula_fixture(), exclude_variables = "L", verbose = FALSE)
    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    fit <- continuous_formula_fit(x, treatment, censor)
    treatment_calls <- treatment$calls()
    expect_gte(length(treatment_calls), 3L)
    rhs <- list(continuous_formula_rhs(treatment_calls[[1L]]$formula),
                continuous_formula_rhs(tail(treatment_calls, 1L)[[1L]]$formula),
                continuous_formula_rhs(continuous_formula_outcome(fit, "event_4")),
                continuous_formula_rhs(continuous_formula_censor(censor, 4L)))
    for (terms in rhs) expect_false(any(grepl("^L_[0-9]+$", terms)))

    # A particular event-indexed value can be excluded without removing the
    # rest of the covariate's history.
    x <- model_formula(continuous_formula_fixture(), exclude_variables = "L_2", verbose = FALSE)
    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    fit <- continuous_formula_fit(x, treatment, censor)
    rhs <- list(continuous_formula_rhs(tail(treatment$calls(), 1L)[[1L]]$formula),
                continuous_formula_rhs(continuous_formula_outcome(fit, "event_4")),
                continuous_formula_rhs(continuous_formula_censor(censor, 4L)))
    for (terms in rhs) {
        expect_false("L_2" %in% terms)
        expect_true("L_3" %in% terms)
    }
})

test_that("marginal censoring applies consistent named and global exclusions", {
    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    x <- model_formula(continuous_formula_fixture(),
                       exclusion_rules = list("^C_[0-9]+$" = "^L_0$"), verbose = FALSE)
    continuous_formula_fit_marginal(x, treatment, censor)
    marginal_rhs <- continuous_formula_rhs(censor$calls()[[1L]]$formula)
    expect_false("L_0" %in% marginal_rhs)
    expect_true("age" %in% marginal_rhs)
    # The censoring-specific rule does not remove L_0 from baseline treatment.
    expect_true("L_0" %in% continuous_formula_rhs(treatment$calls()[[1L]]$formula))

    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    x <- model_formula(continuous_formula_fixture(), exclude_variables = "L", verbose = FALSE)
    continuous_formula_fit_marginal(x, treatment, censor)
    expect_false("L_0" %in% continuous_formula_rhs(censor$calls()[[1L]]$formula))
})

test_that("custom marginal censoring formulas preserve strata and apply rules", {
    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    x <- model_formula(continuous_formula_fixture(), verbose = FALSE)
    continuous_formula_fit_marginal(x, treatment, censor, ~ strata(group) + age)
    default <- censor$calls()[[1L]]$formula
    expect_true("strata(group)" %in% attr(stats::terms(stats::as.formula(default)), "term.labels"))
    expect_true("age" %in% continuous_formula_rhs(default))

    treatment <- continuous_formula_spy()
    censor <- continuous_formula_spy()
    x <- model_formula(continuous_formula_fixture(),
                       exclusion_rules = list("^C_[0-9]+$" = "^age$"),
                       inclusion_rules = list("^C_[0-9]+$" = "L_0"), verbose = FALSE)
    continuous_formula_fit_marginal(x, treatment, censor, ~ strata(group) + age)
    adjusted <- censor$calls()[[1L]]$formula
    expect_true("strata(group)" %in% attr(stats::terms(stats::as.formula(adjusted)), "term.labels"))
    expect_false("age" %in% continuous_formula_rhs(adjusted))
    expect_true("L_0" %in% continuous_formula_rhs(adjusted))

    # Exercise a real Cox fit as well as the spy, including its strata parser.
    x <- model_formula(continuous_formula_fixture(), verbose = FALSE)
    fit <- run_rtmle(x, learner = "learn_glm", estimator = "one-step",
                     censoring_learner = list(name = "cox", fun = "learn_coxph",
                                              marginal_censoring_formula = ~ strata(group) + age))
    expect_true(all(is.finite(fit$estimate$Main_analysis$Estimate)))
})

test_that("one marginal censoring model rejects inconsistent node-specific rules", {
    x <- model_formula(continuous_formula_fixture(),
                       exclusion_rules = list("^C_4$" = "^age$"), verbose = FALSE)
    expect_error(continuous_formula_fit_marginal(x, continuous_formula_spy(),
                                                 continuous_formula_spy()),
                 "[Mm]arginal|shared|different predictors")
})

test_that("custom marginal formulas retain offsets and respect global exclusions", {
    x <- model_formula(continuous_formula_fixture(), verbose = FALSE)
    censor <- continuous_formula_spy()
    continuous_formula_fit_marginal(x, continuous_formula_spy(), censor,
                                    ~ strata(group) + offset(age) + L_0)
    default <- censor$calls()[[1L]]$formula
    expect_match(default, "offset\\(age\\)")
    expect_true("L_0" %in% continuous_formula_rhs(default))
    expect_true("strata(group)" %in% attr(stats::terms(stats::as.formula(default)), "term.labels"))

    x <- model_formula(continuous_formula_fixture(), exclude_variables = "L", verbose = FALSE)
    censor <- continuous_formula_spy()
    continuous_formula_fit_marginal(x, continuous_formula_spy(), censor,
                                    ~ strata(group) + offset(age) + L_0)
    excluded <- censor$calls()[[1L]]$formula
    expect_match(excluded, "offset\\(age\\)")
    expect_false("L_0" %in% continuous_formula_rhs(excluded))
    expect_true("strata(group)" %in% attr(stats::terms(stats::as.formula(excluded)), "term.labels"))

    # Removing a covariate also removes a custom term that depends on its
    # baseline value, instead of leaving strata(L_0) in the fitted model.
    censor <- continuous_formula_spy()
    continuous_formula_fit_marginal(x, continuous_formula_spy(), censor,
                                    ~ strata(L_0) + age)
    excluded_stratum <- censor$calls()[[1L]]$formula
    expect_false("L_0" %in% continuous_formula_rhs(excluded_stratum))
    expect_true("age" %in% continuous_formula_rhs(excluded_stratum))
})

test_that("complex offsets survive unrelated rules and dot keeps formula semantics", {
    x <- model_formula(continuous_formula_fixture(),
                       exclusion_rules = list("^C_[0-9]+$" = "^age$"), verbose = FALSE)
    censor <- continuous_formula_spy()
    continuous_formula_fit_marginal(x, continuous_formula_spy(), censor,
                                    ~ age + offset(log(offset_value)) + strata(group))
    formula <- censor$calls()[[1L]]$formula
    expect_match(formula, "offset\\(log\\(offset_value\\)\\)")
    expect_false("age" %in% continuous_formula_rhs(formula))
    expect_true("strata(group)" %in% attr(stats::terms(stats::as.formula(formula)), "term.labels"))

    # With no applicable rules, preserve the supplied formula so the learner
    # evaluates the dot using its usual formula semantics.
    x <- model_formula(continuous_formula_fixture(), verbose = FALSE)
    censor <- continuous_formula_spy()
    continuous_formula_fit_marginal(x, continuous_formula_spy(), censor, ~ .)
    expect_match(censor$calls()[[1L]]$formula, "~[[:space:]]*\\.")

    # Applicable named rules operate on the expanded formula.
    x <- model_formula(continuous_formula_fixture(),
                       exclusion_rules = list("^C_[0-9]+$" = "^age$"), verbose = FALSE)
    censor <- continuous_formula_spy()
    continuous_formula_fit_marginal(x, continuous_formula_spy(), censor, ~ .)
    rhs <- continuous_formula_rhs(censor$calls()[[1L]]$formula)
    expect_false("age" %in% rhs)
    expect_true("L_0" %in% rhs)
})
