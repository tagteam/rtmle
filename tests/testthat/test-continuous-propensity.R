continuous_propensity_fixture <- function(censored = FALSE) {
    n <- 40L
    a <- rep(c(0, 1, 1, 0), n / 4L)
    d <- data.table::data.table(
        id = seq_len(n), A_0 = a, L_0 = rep(c(0, 1), n / 2L),
        time_0 = 0, event_0 = factor("A"),
        time_1 = rep(c(1, 2), n / 2L),
        event_1 = factor("A"), A_1 = rev(a), L_1 = rep(c(0, 1), n / 2L),
        time_2 = rep(c(3, 4, 5, 6), n / 4L),
        event_2 = factor(if (censored) rep(c("C", "Y", "D", "tauend"), each = n / 4L) else rep(c("Y", "D"), n / 2L)),
        A_2 = rev(a), L_2 = rep(c(0, 1), n / 2L)
    )
    d[, `:=`(time_pooled_2 = time_2, event_pooled_2 = event_2)]
    marginal_data <- d[, .(id, A_0, L_0, time = time_2, event = event_2)]
    list(wide_data = d,
         info = data.table::data.table(time_horizon = 6, last_event = 2L,
                                      is_censored = censored),
         time_covariates = c("A", "L"), baseline_covariates = c("A_0", "L_0"),
         data_marginal_censoring = if (censored) marginal_data else NULL)
}

test_that("continuous propensity reuses learner specifications and leaves inputs alone", {
    p <- continuous_propensity_fixture()
    before <- data.table::copy(p$wide_data)
    seen <- list()
    learner <- function(character_formula, data, intervened_data, value, ...) {
        seen[[length(seen) + 1L]] <<- character_formula
        list(predicted_values = rep(value, nrow(intervened_data)))
    }
    result <- continuous_propensity(p, learner = list(name = "custom", fun = learner, value = 0.3))
    expect_equal(result$data$propensity_0, rep(0.3, 40))
    expect_equal(result$data$propensity_1, rep(0.3, 40))
    expect_length(seen, 2)
    expect_equal(p$wide_data, before)
    expect_equal(result$data$survival_censoring_1, rep(1, 40))
})

test_that("constant treatment is evaluated under the requested intervention", {
    p <- continuous_propensity_fixture()
    p$wide_data[, `:=`(A_0 = 0, A_1 = 0)]
    zero <- continuous_propensity(p, static_intervention = 0)
    one <- continuous_propensity(p, static_intervention = 1)
    expect_equal(zero$data$propensity_0, rep(1, 40))
    expect_equal(zero$data$propensity_1, rep(1, 40))
    expect_equal(one$data$propensity_0, rep(0, 40))
    expect_equal(one$data$propensity_1, rep(0, 40))
})

test_that("propensity bounds apply to newly predicted scores", {
    p <- continuous_propensity_fixture()
    learner <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = rep(0.01, nrow(intervened_data)))
    }
    result <- continuous_propensity(p, learner = learner, gbound = 0.1)
    expect_equal(result$data$propensity_0, rep(0.1, 40))
    expect_equal(result$data$propensity_1, rep(0.1, 40))
    expect_error(continuous_propensity(p, gbound = c(0.1, 0.2)), "gbound")
    expect_error(continuous_propensity(p, gbound = -1), "gbound")
})

test_that("invalid nuisance predictions fail with an informative error", {
    p <- continuous_propensity_fixture()
    bad <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = NA_real_)
    }
    expect_error(continuous_propensity(p, learner = bad), "predictions")
    bad <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = rep(1.1, nrow(intervened_data)))
    }
    expect_error(continuous_propensity(p, learner = bad), "predictions")
    expect_error(continuous_propensity(p, static_intervention = 2), "static_intervention")
    expect_error(continuous_propensity(p, lag = 0), "lag")
    expect_error(continuous_propensity(continuous_propensity_fixture(TRUE), censoring_learner = NULL), "censoring")
})

test_that("Cox interval survival uses the left limit at the endpoint", {
    d <- data.table::data.table(time = c(1, 2, 3, 4), status = c(1, 1, 0, 0))
    pred <- data.table::data.table(time = c(1, 2, 2.5, 3), start = c(0, 0, 1, 2))
    result <- learn_coxph("survival::Surv(time, status) ~ 1", d, pred,
                         time_variable = "time", time_ref = "start")
    expect_equal(result$predicted_values, exp(-c(0, 1/4, 1/3, 0)), tolerance = 1e-12)
    expect_true(all(is.finite(result$predicted_values)))
    expect_error(learn_coxph("survival::Surv(time, status) ~ 1", d,
                            data.table::data.table(time = 1, start = 2),
                            time_ref = "start"), "start")
})

test_that("both censoring strategies return valid interval survival", {
    for (marginal in c(TRUE, FALSE)) {
        p <- continuous_propensity_fixture(censored = TRUE)
        censor <- list(name = "coxph", fun = "learn_coxph", marginal_censoring = marginal)
        result <- continuous_propensity(p, censoring_learner = censor)
        expect_true(all(result$data$survival_censoring_2 > 0))
        expect_true(all(result$data$survival_censoring_2 <= 1))
        expect_equal(result$data$survival_censoring_2,
                     result$data$survival_censoring_pooled_2)
    }
})

test_that("continuous nuisance fits reproduce the source seed-34 oracle", {
    # These values come from contICEIPCW's existing n = 1000 conservative
    # regression tests. The source package is not needed when running them.
    set.seed(34)
    d <- simulate_continuous_time_data(n = 1000, no_competing_events = TRUE,
                                       uncensored = FALSE)
    expected <- list(
        marginal = c(0.986264290518851, 0.984672060697723, 0.972597415471961),
        conditional = c(0.986155224477072, 0.987601297850103, 0.969323799924238)
    )
    for (marginal in c(TRUE, FALSE)) {
        p <- prepare_continuous_data(d, time_horizons = 720,
                                      time_covariates = c("A", "L"),
                                      baseline_covariates = c("age", "A_0", "L_0"))
        censor <- list(name = "coxph", fun = "learn_coxph", marginal_censoring = marginal)
        result <- continuous_propensity(p, censoring_learner = censor)
        survival_means <- vapply(paste0("survival_censoring_", 1:3),
                                 function(column) mean(result$data[[column]], na.rm = TRUE),
                                 numeric(1))
        expect_equal(unname(survival_means),
                     expected[[if (marginal) "marginal" else "conditional"]],
                     tolerance = 1e-10)
        expect_equal(mean(result$data$propensity_1^2, na.rm = TRUE),
                     0.360399927683227, tolerance = 1e-10)
        expect_equal(mean(result$data$propensity_2^2, na.rm = TRUE),
                     0.312896735907919, tolerance = 1e-10)
    }
})

test_that("censoring learner warnings are retained in diagnostics", {
    p <- continuous_propensity_fixture(censored = TRUE)
    warning_learner <- function(character_formula, data, intervened_data, ...) {
        warning("censoring model warning for review")
        list(predicted_values = rep(0.9, nrow(intervened_data)))
    }
    result <- NULL
    expect_warning(result <- continuous_propensity(p, censoring_learner = list(name = "warning", fun = warning_learner, marginal_censoring = FALSE)), NA)
    expect_true(any(grepl("censoring model warning for review",
                           result$diagnostics$learner_warnings$warning, fixed = TRUE)))
})
