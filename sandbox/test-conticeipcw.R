# Standalone testthat comparisons of RTMLE and contICEIPCW.
# Run from the RTMLE source directory with contICEIPCW installed:
# Rscript --vanilla -e 'pkgload::load_all(".", quiet = TRUE, export_all = FALSE); testthat::test_file("sandbox/test-conticeipcw.R", stop_on_failure = TRUE)' </dev/null

testthat::local_edition(3)

# IPC weights and treatment probabilities

parity_weights_fixture <- function(marginal, horizon = 720, treatment = 1,
                                   lag = NULL, exclude_latest_covariate = NULL) {
    testthat::skip_if_not_installed("contICEIPCW")
    data <- withr::with_seed(3901, rtmle::simulate_continuous_time_data(
        750, static_intervention_baseline = NULL,
        baseline_rate_list = list(A = .005, L = .001, C = .0004, Y = .0001, D = .00015)
    ))
    time_covariates <- c("A", "L")
    baseline_covariates <- c("age", "A_0", "L_0")
    reference <- contICEIPCW::prepare_data(
        lapply(data, data.table::copy), horizon, time_covariates, baseline_covariates,
        marginal_censoring = marginal
    )
    target <- rtmle:::prepare_continuous_data(data, horizon, time_covariates, baseline_covariates)
    reference <- contICEIPCW::propensity_scores(
        reference, model_treatment = "learn_glm_logistic", model_hazard = "learn_coxph",
        static_intervention = treatment, lag = lag,
        exclude_latest_covariate = exclude_latest_covariate
    )
    target <- rtmle:::continuous_propensity(
        target, learner = "learn_glm", static_intervention = treatment, lag = lag,
        exclude_latest_covariate = exclude_latest_covariate, keep_censoring_prediction_models = TRUE,
        censoring_learner = list(name = "coxph", fun = "learn_coxph", marginal_censoring = marginal)
    )
    testthat::expect_identical(target$data$id, reference$data$id)
    last_event <- target$prepared_data_object$info$last_event
    testthat::expect_equal(last_event, 3)
    testthat::expect_equal(last_event, reference$prepared_data_object$info$last_event)

    for (k in seq_len(last_event) - 1L) {
        column <- paste0("propensity_", k)
        observed <- target$data[[column]]
        decisions <- target$data[[paste0("event_", k)]] %in% "A"
        info <- sprintf("Treatment %d at event index %d", treatment, k)
        testthat::expect_true(sum(decisions) > 1L, info = info)
        testthat::expect_true(stats::sd(observed[decisions]) > 0, info = info)
        testthat::expect_identical(!is.na(observed), decisions, info = info)
        testthat::expect_equal(observed, reference$data[[column]], tolerance = 1e-10, info = info)
    }
    list(reference = reference, target = target, horizon = horizon)
}

parity_weights_compare_ipcw <- function(fixture) {
    target <- fixture$target
    reference <- fixture$reference
    horizon <- fixture$horizon
    last_event <- target$prepared_data_object$info$last_event
    fit <- rtmle:::continuous_sequential_regression(
        target, th = horizon, a = 1,
        learner = rtmle:::continuous_nuisance_learner("learn_ipcw_glm"),
        estimator = "sequential_regression", lag = NULL, seed = NULL,
        prediction_range = NULL, stratify = FALSE, intervene_all = FALSE,
        intervention_hook = NULL, save_fitted_objects = TRUE
    )

    for (k in seq_len(last_event)) {
        previous_event <- paste0("event_", k - 1L)
        previous_time <- paste0("time_", k - 1L)
        survival <- paste0("survival_censoring_", k)
        testthat::expect_equal(target$data[[survival]], reference$data[[survival]],
                     tolerance = 1e-11, info = paste("Interval survival at event index", k))
        at_risk <- target$data[[previous_event]] %in% c("A", "L") &
            target$data[[previous_time]] < horizon
        reference_rows <- reference$data[which(at_risk)]
        data.table::set(reference_rows, j = "event_k_prev", value = reference_rows[[previous_event]])
        if (k == last_event) {
            for (variable in c("time", "event", "survival_censoring")) {
                data.table::set(reference_rows, j = paste0(variable, "_", k),
                                value = reference_rows[[paste0(variable, "_pooled_", k)]])
            }
        }
        expected <- contICEIPCW:::ipcw_k(
            reference_rows, k, reference$marginal_censoring_fit, horizon,
            is_censored = TRUE, survival_function = reference_rows[[survival]]
        )
        actual <- stats::weights(fit$models[[paste0("event_", k)]]$fit, type = "prior")
        testthat::expect_equal(unname(actual), expected, tolerance = 1e-11,
                     info = paste("Regression IPCW at event index", k))
        testthat::expect_true(any(actual == 0))
        testthat::expect_true(any(actual > 1))
    }
    invisible(fit)
}

testthat::test_that("event-indexed Cox IPCW agrees with contICEIPCW through end of follow-up", {
    # The source's event-indexed IPCW uses its fast path, which requires the
    # horizon to cover all event times for comparison of all regression weights.
    fixture <- parity_weights_fixture(marginal = FALSE, horizon = 900)
    testthat::expect_true(all(fixture$target$data$time_3 <= fixture$horizon, na.rm = TRUE))
    parity_weights_compare_ipcw(fixture)
})

testthat::test_that("marginal Cox IPCW agrees with contICEIPCW at events and at the horizon", {
    fixture <- parity_weights_fixture(marginal = TRUE)
    fit <- parity_weights_compare_ipcw(fixture)
    data <- fixture$target$data
    at_risk <- which(data$event_2 %in% c("A", "L") & data$time_2 < fixture$horizon)
    later <- data$time_pooled_3[at_risk] > fixture$horizon
    testthat::expect_true(any(later))
    weights <- stats::weights(fit$models$event_3$fit, type = "prior")
    testthat::expect_true(any(weights[later] > 1))
})

testthat::test_that("logistic propensity scores agree with contICEIPCW for each static treatment", {
    for (treatment in 0:1) {
        parity_weights_fixture(marginal = TRUE, treatment = treatment,
                               lag = 1L, exclude_latest_covariate = "time")
    }
})

testthat::test_that("cumulative treatment and censoring weights agree with contICEIPCW on the risk set", {
    for (marginal in c(FALSE, TRUE)) {
        fixture <- parity_weights_fixture(marginal = marginal)
        target <- fixture$target$data
        reference <- data.table::copy(fixture$reference$data)
        last_event <- fixture$target$prepared_data_object$info$last_event
        # The reference cumulative helper multiplies treatment values directly,
        # so its cumulative adherence is valid for the always-treated regime.
        contICEIPCW:::cumulative_inverse_probability_weights(
            reference, static_intervention = 1, time_horizon = fixture$horizon,
            return_ipw = FALSE, last_event = last_event
        )
        actual <- rtmle:::continuous_weights(target, a = 1, last_event = last_event)
        for (k in seq_len(last_event)) {
            at_risk <- target[[paste0("event_", k - 1L)]] %in% c("A", "L")
            expected <- reference[[paste0("ipw_cum_weight_", k - 1L)]]
            testthat::expect_equal(actual[at_risk, k], expected[at_risk], tolerance = 1e-10,
                         info = sprintf("Marginal Cox: %s, event index %d", marginal, k))
            testthat::expect_true(any(actual[at_risk, k] > 0))
            testthat::expect_true(any(actual[at_risk, k] == 0))
            testthat::expect_true(all(actual[!at_risk, k] == 0))
        }
    }
})

# Backward outcome regression and the compiled solver

parity_regression_capture <- function(outcome) {
    stages <- list()
    list(
        hook = function(data, k) {
            if (outcome %in% names(data)) {
                rows <- order(data$id)
                stages[[paste0("event_", k + 1L)]] <<- data.frame(
                    id = data$id[rows],
                    treatment = as.numeric(as.character(data[[paste0("A_", k)]][rows])),
                    pseudo_outcome = data[[outcome]][rows]
                )
            }
            data
        },
        stages = function() stages
    )
}

testthat::test_that("continuous-time backward regressions agree at each event index", {
    testthat::skip_if_not_installed("contICEIPCW")
    withr::local_seed(34)
    horizon <- 720
    data <- rtmle::simulate_continuous_time_data(600, static_intervention_baseline = NULL)
    preparation <- list(time_horizons = horizon, time_covariates = c("A", "L"),
                        baseline_covariates = c("age", "A_0", "L_0"))
    prepared <- do.call(rtmle:::prepare_continuous_data, c(list(data = data), preparation))
    reference_prepared <- do.call(contICEIPCW::prepare_data,
        c(list(data = data.table::copy(data), marginal_censoring = TRUE), preparation))
    reference <- contICEIPCW::propensity_scores(
        reference_prepared, model_treatment = "learn_glm_logistic", model_hazard = "learn_coxph"
    )
    nuisance <- rtmle:::continuous_propensity(
        prepared, learner = NULL, keep_censoring_prediction_models = TRUE
    )

    # Capture the reference predictions as the backward recursion uses them.
    reference_fit <- contICEIPCW::debias_ice_ipcw
    reference_predict <- contICEIPCW:::predict_intervention
    environment(reference_fit) <- list2env(list(
        predict_intervention = function(data, k, ...) {
            values <- reference_predict(data, k, ...)
            # The reference also predicts at baseline again for the final mean.
            if ("pseudo_outcome" %in% names(data)) {
                rows <- order(data$id)
                expected_predictions[[paste0("event_", k + 1L)]] <<-
                    stats::setNames(as.numeric(values)[rows], data$id[rows])
            }
            values
        }
    ), parent = environment(reference_fit))

    learners <- c(oipcw_expit = "learn_oipcw", ipcw_glm_expit = "learn_ipcw_glm")
    testthat::expect_equal(prepared$info$last_event, 3)
    for (learner in names(learners)) {
        observed_stages <- parity_regression_capture("rtmle_pseudo_outcome")
        expected_stages <- parity_regression_capture("pseudo_outcome")
        expected_predictions <- list()
        # A grid with only its endpoints retains the horizon predictions.
        fit <- rtmle:::continuous_sequential_regression(
            nuisance = nuisance, th = horizon, a = 1,
            learner = rtmle:::continuous_nuisance_learner(learners[[learner]]),
            estimator = "sequential_regression", lag = NULL, seed = NULL,
            prediction_range = c(0, 1), stratify = FALSE, intervene_all = FALSE,
            intervention_hook = observed_stages$hook, save_fitted_objects = TRUE,
            q_grid = c(0, horizon)
        )
        # contICEIPCW obtains the plug-in risk within its one-step workflow.
        expected <- reference_fit(
            reference, model_pseudo_outcome = learner, conservative = TRUE,
            intervention_hook = expected_stages$hook, save_coefficients = TRUE
        )
        stages <- observed_stages$stages()
        testthat::expect_equal(stages, expected_stages$stages(),
                               tolerance = 1e-7, info = learner)
        expected_coefficients <- expected$coefficients_pseudo_outcome[[paste0("time_horizon_", horizon)]]
        testthat::expect_length(fit$models, length(expected_coefficients))
        for (k in seq_along(expected_coefficients)) {
            event <- paste0("event_", k)
            ids <- stages[[event]]$id
            rows <- match(ids, nuisance$data$id)
            predictions <- stats::setNames(fit$q_grid[[k]][rows, 2L], ids)
            testthat::expect_equal(predictions, expected_predictions[[event]],
                                   tolerance = 1e-7, info = paste(learner, event))

            summary <- fit$models[[event]]$fit_summary
            coefficients <- stats::setNames(summary[, 1], rownames(summary))
            # QR aliases are zero-filled by contICEIPCW and NA in a fitted GLM.
            coefficients[is.na(coefficients)] <- 0
            comparison <- expected_coefficients[[k]]
            unused <- setdiff(names(coefficients), names(comparison))
            testthat::expect_equal(unname(coefficients[unused]), numeric(length(unused)))
            testthat::expect_equal(coefficients[names(comparison)], comparison,
                                   tolerance = 1e-7, info = paste(learner, "event", k))
        }
        testthat::expect_equal(fit$estimate, expected$results$ice_ipcw_estimate,
                               tolerance = 1e-7, info = learner)
    }
})

testthat::test_that("the C++ estimating-equation solvers agree on IPCW pseudo-outcomes", {
    testthat::skip_if_not_installed("contICEIPCW")
    withr::local_seed(314)
    n <- 300
    history <- cbind(1, stats::rnorm(n), stats::rbinom(n, 1, 0.5))
    offset <- stats::rnorm(n, sd = 0.2)
    risk <- stats::plogis(history %*% c(-1.5, 0.35, -0.2) + offset)
    event <- stats::rbinom(n, 1, risk)
    survival <- stats::runif(n, 0.67, 0.97)
    observed <- stats::rbinom(n, 1, survival)
    outcome <- observed * event / survival
    beta <- c(-1, 0.1, 0.2)
    testthat::expect_true(any(outcome > 1))

    for (weights in list(NULL, observed / survival)) {
        fit <- rtmle:::estimating_equation_cpp(history, outcome, beta, offset, weights)
        expected <- contICEIPCW:::estimating_equation_cpp(history, outcome, beta, offset, weights)
        testthat::expect_equal(as.numeric(fit), as.numeric(expected), tolerance = 1e-8)

        residual <- outcome - stats::plogis(history %*% fit + offset)
        if (!is.null(weights)) residual <- weights * residual
        testthat::expect_lt(max(abs(crossprod(history, residual))), 1e-7)
    }
})

# Individual nuisance learners
parity_learners_data <- function() {
    withr::with_seed(271, {
        n <- 400L
        L <- stats::rnorm(n)
        Z <- stats::rnorm(n)
        A <- stats::rbinom(n, 1, stats::plogis(0.4 * L - 0.3 * Z))
        G <- stats::plogis(1.3 + 0.2 * L + 0.25 * A)
        observed <- stats::rbinom(n, 1, G)
        risk <- stats::plogis(-1.1 + 0.3 * A - 0.3 * L)
        event <- stats::rbinom(n, 1, risk)
        continuation <- stats::rbinom(n, 1, 0.4)
        outcome <- ifelse(continuation == 1, risk, event)
        ipcw <- observed / G
        censoring_time <- stats::rexp(n, rate = exp(-0.3 + 0.2 * L + 0.3 * A))
        failure_time <- stats::rexp(n, rate = exp(0.25 + 0.4 * Z))
        data.table::data.table(id = seq_len(n), A = A, L = L, Z = Z,
                               outcome = outcome, weight = ipcw * outcome, ipcw = ipcw,
                               C = as.integer(censoring_time <= failure_time),
                               time = pmin(censoring_time, failure_time))
    })
}

parity_learners_reference_q <- function(learner, data) {
    fit <- getFromNamespace("learn_Q", "contICEIPCW")
    fit(model_type = learner, history_of_variables = c("A", "L"),
        data_learn = data.table::copy(data), outcome_name = "weight",
        outcome_string_unweighted = "outcome", ipcw_name = "ipcw",
        penalize = FALSE, verbose = FALSE, k = 1L, save_coefficients = TRUE)
}

testthat::test_that("logistic treatment learners give the same fitted propensities", {
    testthat::skip_if_not_installed("contICEIPCW")
    d <- parity_learners_data()
    formula <- "A ~ L + Z"
    fit_reference <- getFromNamespace("learn_glm_logistic", "contICEIPCW")

    reference <- fit_reference(formula, d, penalize = FALSE, save_coefficients = TRUE)
    actual <- rtmle::learn_glm(formula, d, d, save_fitted_objects = TRUE)
    testthat::expect_equal(as.numeric(actual$predicted_values), reference$predictions, tolerance = 1e-7)
    testthat::expect_equal(stats::coef(actual$fit), reference$coefficients, tolerance = 1e-6)

    # Both fits select lambda.min on the same folds. RTMLE supplies an
    # intercept column that glmnet excludes, so allow for solver precision.
    reference <- withr::with_seed(567, fit_reference(formula, d, penalize = TRUE))
    actual <- withr::with_seed(567, rtmle::learn_glmnet(formula, d, d, selector = "min",
                                              family = "binomial"))
    testthat::expect_equal(actual$predicted_values, reference$predictions, tolerance = 1e-5)
})

testthat::test_that("Cox censoring learners give the same left-limit survival predictions", {
    testthat::skip_if_not_installed("contICEIPCW")
    d <- parity_learners_data()
    formula <- "survival::Surv(time, C) ~ A + L"
    fit_reference <- getFromNamespace("learn_coxph", "contICEIPCW")
    reference <- fit_reference(formula, d, penalize = FALSE)
    actual <- rtmle::learn_coxph(formula, d, left_limit = TRUE)

    testthat::expect_true(any(reference$pred < 1))
    testthat::expect_equal(actual$predicted_values, reference$pred, tolerance = 1e-10)
    testthat::expect_equal(stats::coef(actual$fit), stats::coef(reference$fit), tolerance = 1e-10)
})

testthat::test_that("IPCW GLM learners agree for logit and probit links", {
    testthat::skip_if_not_installed("contICEIPCW")
    d <- parity_learners_data()
    intervened <- d[30:1, .(A = 1L, L)]
    reference_logit <- parity_learners_reference_q("ipcw_glm_expit", d)
    # The reference binomial probit fit warns about fractional weighted responses.
    testthat::expect_warning(
        reference_probit <- parity_learners_reference_q("ipcw_glm_probit", d),
        "non-integer #successes"
    )

    for (link in c("logit", "probit")) {
        reference <- if (link == "logit") reference_logit else reference_probit
        actual <- rtmle::learn_ipcw_glm("outcome ~ A + L", d, intervened, weights = d$ipcw,
                                 link = link, save_fitted_objects = TRUE)
        testthat::expect_equal(actual$predicted_values, as.numeric(reference$predict_fun(intervened)),
                     tolerance = 1e-10, info = link)
        testthat::expect_equal(stats::coef(actual$fit), reference$coefficients,
                     tolerance = 1e-10, info = link)
    }
})

testthat::test_that("bounded OIPCW learners agree for responses above one", {
    testthat::skip_if_not_installed("contICEIPCW")
    d <- parity_learners_data()
    intervened <- d[30:1, .(A = 1L, L)]
    testthat::expect_true(any(d$weight > 1))
    reference <- parity_learners_reference_q("oipcw_expit", d)
    actual <- rtmle::learn_oipcw("weight ~ A + L", d, intervened, save_fitted_objects = TRUE)

    testthat::expect_equal(actual$predicted_values, as.numeric(reference$predict_fun(intervened)),
                 tolerance = 1e-8)
    testthat::expect_equal(actual$fit$coefficients, reference$coefficients, tolerance = 1e-8)
})

testthat::test_that("the linear outcome model agrees with a Gaussian learner at lambda zero", {
    testthat::skip_if_not_installed("contICEIPCW")
    d <- parity_learners_data()
    intervened <- d[30:1, .(A = 0L, L)]
    reference <- parity_learners_reference_q("lm", d)
    actual <- rtmle::learn_glmnet("weight ~ A + L", d, intervened, family = "gaussian",
                           lambda = 0)
    testthat::expect_equal(actual$predicted_values, as.numeric(reference$predict_fun(intervened)),
                 tolerance = 1e-7)
})

# One-step estimation

parity_one_step_data <- function(uncensored) {
    withr::with_seed(34L, rtmle::simulate_continuous_time_data(
        n = 900L, uncensored = uncensored,
        static_intervention_baseline = NULL,
        baseline_rate_list = list(A = 0.005, L = 0.001, C = 0.0004, Y = 0.0001, D = 0.00015)
    ))
}

parity_one_step_analysis <- function(data, horizon, treatment) {
    events <- data$timevarying_data
    event_rows <- function(marks, variables) {
        events[events$event %in% marks, variables, with = FALSE]
    }
    updates <- function(mark) {
        out <- event_rows(mark, c("id", "time", mark))
        data.table::setnames(out, mark, "value")
        out
    }
    x <- rtmle::rtmle_init(c(0, horizon), name_id = "id", name_time = "time",
                          name_outcome = "Y", name_competing = "D",
                          name_censoring = "C", continuous_time = TRUE)
    x <- rtmle::add_baseline_data(x, data$baseline_data)
    x <- rtmle::add_continuous_data(
        x, outcome_data = event_rows("Y", c("id", "time")),
        competing_data = event_rows("D", c("id", "time")),
        censored_data = event_rows(c("C", "tauend"), c("id", "time", "event")),
        timevar_data = list(A = updates("A"), L = updates("L")),
        treatment_events = "A", covariate_events = "L"
    )
    x <- rtmle::prepare_rtmle_data(x)
    x <- rtmle::regime(x, "static", intervention = treatment, treatment_variables = "A")
    x <- rtmle::target(x, "risk", regimes = "static")
    rtmle::model_formula(x, verbose = FALSE)
}

parity_one_step_fit <- function(data, treatment, marginal = TRUE) {
    horizon <- 720
    reference_data <- data.table::copy(data)
    if (treatment == 0) {
        # The reference cumulative weight helper multiplies treatment values,
        # so encode the desired arm as 1 throughout the reference data.
        reference_data$baseline_data$A_0 <- 1 - reference_data$baseline_data$A_0
        reference_data$timevarying_data$A <- 1 - reference_data$timevarying_data$A
    }
    reference <- contICEIPCW::prepare_data(
        reference_data, time_horizons = horizon,
        time_covariates = c("A", "L"), baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = marginal
    )
    reference <- contICEIPCW::propensity_scores(
        reference, model_treatment = "learn_glm_logistic", model_hazard = "learn_coxph",
        static_intervention = 1
    )
    ref <- contICEIPCW::debias_ice_ipcw(
        reference, model_pseudo_outcome = "oipcw_expit", model_hazard = "learn_coxph",
        conservative = TRUE, static_intervention = 1, return_ic = TRUE
    )

    x <- parity_one_step_analysis(data, horizon, treatment)
    # The working checkout contains a browser() call. In a non-interactive
    # Rscript process with closed stdin, allow R to continue past it.
    withr::local_envvar(c("_R_CHECK_BROWSER_NONINTERACTIVE_" = "false"))
    invisible(utils::capture.output(fit <- rtmle::run_rtmle(
        x, estimator = "one-step", time_horizon = 1,
        learner = "learn_oipcw", propensity_learner = "learn_glm",
        censoring_learner = list(name = "coxph", fun = "learn_coxph",
                                 marginal_censoring = marginal),
        conservative = TRUE
    )))
    fit_ids <- fit$continuous_nuisance$static$data$id
    reference_order <- match(fit_ids, reference$data$id)
    testthat::expect_setequal(fit_ids, reference$data$id)
    actual <- fit$estimate$Main_analysis
    expected <- ref$results
    list(
        actual = list(estimate = actual$Estimate, se = actual$Standard_error,
                      ice = actual$ICE_IPCW, ipw = actual$IPW,
                      ic = fit$IC$risk$static$time_horizon_1),
        expected = list(estimate = expected$estimate, se = expected$se,
                        ice = expected$ice_ipcw_estimate, ipw = expected$ipw,
                        ic = ref$ic$ic[[1L]][reference_order])
    )
}

testthat::test_that("conservative one-step updates agree for both static regimes", {
    testthat::skip_if_not_installed("contICEIPCW")
    uncensored <- parity_one_step_data(TRUE)
    censored <- parity_one_step_data(FALSE)
    testthat::expect_true(any(censored$timevarying_data$event == "C" &
                    censored$timevarying_data$time <= 720))
    testthat::expect_true(any(censored$timevarying_data$event == "D" &
                    censored$timevarying_data$time <= 720))

    for (treatment in 0:1) {
        fit <- parity_one_step_fit(uncensored, treatment)
        testthat::expect_equal(fit$actual, fit$expected, tolerance = 1e-8, info = paste("No dropout,", treatment))
        for (marginal in c(FALSE, TRUE)) {
            fit <- parity_one_step_fit(censored, treatment, marginal)
            testthat::expect_equal(fit$actual, fit$expected, tolerance = 1e-8,
                         info = sprintf("Marginal Cox: %s, treatment: %d", marginal, treatment))
            adjustment <- mean(fit$actual$ic)
            testthat::expect_gt(abs(adjustment), 1e-5)
            testthat::expect_equal(fit$actual$estimate, fit$actual$ice + adjustment, tolerance = 1e-8)
        }
    }
})
