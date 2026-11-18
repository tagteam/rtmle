test_that("stratified Cox survival uses each stratum's hazard and endpoint left limit", {
    d <- data.frame(time = rep(1:4, 2), status = c(1, 1, 0, 0, 1, 0, 1, 0),
                    group = rep(c("a", "b"), each = 4))
    p <- data.frame(time = c(1, 2, 3, 4), start = c(0, 0, 1, 1),
                    group = c("a", "a", "b", "b"))
    fit <- learn_coxph("survival::Surv(time, status) ~ strata(group)", d, p,
                       time_variable = "time", time_ref = "start")
    expect_equal(fit$predicted_values, exp(-c(0, 1/4, 0, 1/2)), tolerance = 1e-12)
    expect_error(learn_coxph("survival::Surv(time, status) ~ strata(group)", d,
        data.frame(time = 2, start = 0, group = "unknown"), time_ref = "start"), "level|strat")
})

test_that("censoring learners reject left-censored survival responses", {
    d <- data.frame(time = 1:4, status = c(1, 0, 1, 0), group = c(0, 1, 0, 1))
    formula <- "survival::Surv(time, status, type = 'left') ~ group"
    for (learner in list(learn_coxph, learn_weibull)) {
        expect_error(learner(formula, d), "right-censored")
        d$status <- rep(0, nrow(d))
        expect_error(learner(formula, d), "right-censored")
        d$status <- c(1, 0, 1, 0)
    }
})

test_that("Weibull predictions equal the fitted survival ratio and reuse the fit", {
    d <- data.frame(time = c(1, 2, 3, 4, 5, 6, 7, 8), status = c(1, 1, 0, 1, 0, 1, 0, 1),
                    group = rep(c(0, 1), 4))
    p <- data.frame(time = c(1, 4, 7), start = c(0, 2, 7), group = c(0, 1, 0))
    fit <- learn_weibull("survival::Surv(time, status) ~ group", d, p, time_ref = "start")
    direct <- survival::survreg(survival::Surv(time, status) ~ group, d, dist = "weibull")
    location <- as.numeric(stats::predict(direct, p, type = "lp"))
    s <- function(t) exp(-(t / exp(location))^(1 / direct$scale))
    expect_equal(fit$predicted_values, s(p$time) / s(p$start), tolerance = 1e-12)
    again <- learn_weibull("survival::Surv(time, status) ~ group", d, p,
                            time_ref = "start", reuse_fit = fit)
    expect_equal(again$predicted_values, fit$predicted_values)
    expect_equal(fit$predicted_values[3], 1)
    expect_error(learn_weibull("survival::Surv(time, status) ~ group", d,
                              data.frame(time = 1, start = 2, group = 0), time_ref = "start"), "start")
    d$status <- 0
    expect_equal(learn_weibull("survival::Surv(time, status) ~ group", d, p,
                               time_ref = "start")$predicted_values, rep(1, nrow(p)))
})

censoring_learner_fixture <- function() {
    n <- 8L
    d <- data.table::data.table(id = 1:n, A_0 = rep(c(0, 1), 4), L_0 = rep(c(0, 1), each = 4),
        time_0 = 0, event_0 = "A", time_1 = rep(c(1, 2), 4), event_1 = "A",
        A_1 = rep(c(0, 1), 4), L_1 = rep(c(0, 1), each = 4),
        time_2 = rep(c(3, 4, 5, 6), 2), event_2 = rep(c("C", "Y", "C", "Y"), 2))
    d[, `:=`(time_pooled_2 = time_2, event_pooled_2 = event_2)]
    list(wide_data = d, info = data.table::data.table(time_horizon = 6, last_event = 2L, is_censored = TRUE),
         time_covariates = c("A", "L"), baseline_covariates = c("A_0", "L_0"),
         data_marginal_censoring = d[, .(id, A_0, L_0, time = time_2, event = event_2)])
}

test_that("marginal censoring can reuse an arbitrary survival learner", {
    fit_calls <- 0L
    reuse_calls <- 0L
    exponential <- function(character_formula, data, intervened_data, time_variable,
                            time_ref = NULL, reuse_fit = NULL, ...) {
        if (is.null(reuse_fit)) fit_calls <<- fit_calls + 1L else reuse_calls <<- reuse_calls + 1L
        start <- if (is.null(time_ref)) 0 else intervened_data[[time_ref]]
        list(predicted_values = exp(-0.1 * (intervened_data[[time_variable]] - start)), rate = 0.1)
    }
    treatment <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = rep(0.5, nrow(intervened_data)))
    }
    p <- censoring_learner_fixture()
    result <- continuous_propensity(p, learner = treatment, censoring_learner = exponential)
    expect_equal(fit_calls, 1L)
    expect_gt(reuse_calls, 0L)
    expect_equal(result$data$survival_censoring_1, exp(-0.1 * p$wide_data$time_1))
    expect_equal(result$data$survival_censoring_2, exp(-0.1 * (p$wide_data$time_2 - p$wide_data$time_1)))
    result2 <- continuous_propensity(p, learner = treatment, censoring_learner = "learn_weibull")
    expect_true(all(is.finite(result2$data$survival_censoring_1)))
    result3 <- continuous_propensity(p, learner = treatment, censoring_learner = list(name = "coxph", fun = "learn_coxph", marginal_censoring_formula = ~ strata(L_0)))
    expect_true(all(is.finite(result3$data$survival_censoring_1)))
})

test_that("marginal censoring is specified with the censoring learner", {
    marginal <- list(name = "coxph", fun = "learn_coxph", marginal_censoring = TRUE)
    conditional <- list(name = "coxph", fun = "learn_coxph", marginal_censoring = FALSE)
    p <- censoring_learner_fixture()
    first <- continuous_propensity(p, censoring_learner = marginal)
    second <- continuous_propensity(p, censoring_learner = conditional)
    expect_false(is.null(first$models$censoring$marginal))
    expect_true(is.null(second$models$censoring$marginal))
    expect_true(all(is.finite(first$data$survival_censoring_1)))
    expect_true(all(is.finite(second$data$survival_censoring_1)))
    expect_error(continuous_propensity(p, censoring_learner = list(
        name = "coxph", fun = "learn_coxph", marginal_censoring = "yes")), "marginal_censoring")
    expect_error(continuous_propensity(p, censoring_learner = list(
        name = "coxph", fun = "learn_coxph", marginal_censoring = FALSE,
        marginal_censoring_formula = ~ L_0)), "requires marginal_censoring")
    expect_error(continuous_propensity(p, censoring_learner = list(
        name = "coxph", fun = "learn_coxph", marginal_censoring_formula = "not a formula")), "one-sided")
})

test_that("one prepared analysis can switch censoring model structures", {
    set.seed(34)
    d <- simulate_continuous_time_data(300, uncensored = FALSE)
    x <- continuous_analysis(d, 720, regimes = "Always_A", min_events = 10)
    marginal <- run_rtmle(x, learner = "learn_oipcw", propensity_learner = "glm",
                          estimator = "one-step")
    conditional_learner <- list(name = "coxph", fun = "learn_coxph", marginal_censoring = FALSE)
    switched <- run_rtmle(marginal, learner = "learn_oipcw", propensity_learner = "glm",
                           censoring_learner = conditional_learner, estimator = "one-step", refit = FALSE)
    fresh <- run_rtmle(x, learner = "learn_oipcw", propensity_learner = "glm",
                       censoring_learner = conditional_learner, estimator = "one-step")
    expect_equal(switched$estimate$Main_analysis, fresh$estimate$Main_analysis)
    expect_false(isTRUE(all.equal(marginal$continuous_nuisance$Always_A$data$survival_censoring_1,
                                 switched$continuous_nuisance$Always_A$data$survival_censoring_1)))
    expect_error(prepare_rtmle_data(x, marginal_censoring = FALSE), "marginal_censoring")
})
