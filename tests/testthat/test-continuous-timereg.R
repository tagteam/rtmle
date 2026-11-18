test_that("Aalen intercept agrees with Nelson-Aalen interval survival", {
    skip_if_not_installed("timereg")
    d <- data.frame(id = 1:6, time = 1:6, status = c(1, 0, 1, 0, 1, 0))
    p <- data.frame(time = c(1, 3, 5, 6), start = 0)
    fit <- learn_aalen("survival::Surv(time, status) ~ 1", d, p, time_ref = "start")
    expect_equal(fit$predicted_values, exp(-c(0, 1/6, 1/6 + 1/4, 1/6 + 1/4 + 1/2)))
    p$start <- c(0, 1, 3, 5)
    reused <- learn_aalen("ignored", d, p, time_ref = "start", reuse_fit = fit)
    expect_equal(reused$predicted_values, c(1, 1, 1, 1))
    right <- learn_aalen("ignored", d, p, time_ref = "start", reuse_fit = fit, left_limit = FALSE)
    expect_equal(right$predicted_values, exp(-c(1/6, 1/4, 1/2, 0)))
})

test_that("time-varying Aalen strata agree with manually counted risk sets", {
    skip_if_not_installed("timereg")
    n <- 8
    d <- data.frame(id = rep(1:n, each = 2), start = rep(c(0, 1), n),
                    time = as.vector(rbind(rep(1, n), 2 + seq_len(n))),
                    status = as.vector(rbind(rep(0, n), c(rep(1, 4), rep(0, 4)))),
                    group = factor(as.vector(rbind(rep(c("a", "b"), 4), rep(c("b", "a"), 4)))))
    p <- data.frame(group = factor(c("a", "b"), levels = c("a", "b")), start = 1, time = 7)
    fit <- learn_aalen("survival::Surv(start, time, status) ~ -1 + group", d, p, time_ref = "start")
    expected <- vapply(as.character(p$group), function(group) {
        hazard <- 0
        for (u in d$time[d$status == 1]) {
            risk <- d$start < u & d$time >= u & d$group == group
            fails <- d$time == u & d$status == 1 & d$group == group
            hazard <- hazard + sum(fails) / sum(risk)
        }
        exp(-hazard)
    }, numeric(1))
    expect_equal(fit$predicted_values, expected, ignore_attr = TRUE)
    single <- learn_aalen("ignored", d, p[1, , drop = FALSE], time_ref = "start", reuse_fit = fit)
    expect_equal(single$predicted_values, unname(expected[1]))
})

test_that("Cox-Aalen proportional model agrees with independently fitted Cox", {
    skip_if_not_installed("timereg")
    set.seed(534)
    d <- data.frame(id = 1:100, time = rexp(100), status = rbinom(100, 1, .7), A = rbinom(100, 1, .5))
    p <- data.frame(time = c(.5, 1, 1.5), start = c(0, .2, .5), A = c(0, 1, 0))
    fitted <- learn_cox_aalen("survival::Surv(time, status) ~ prop(A)", d, p, time_ref = "start")
    oracle <- survival::coxph(survival::Surv(time, status) ~ A, data = d, ties = "breslow")
    baseline <- survival::basehaz(oracle, centered = FALSE)
    h <- function(t) c(0, baseline$hazard)[findInterval(t, baseline$time) + 1L]
    expected <- exp(-exp(p$A * stats::coef(oracle)) * (h(p$time) - h(p$start)))
    expect_equal(fitted$predicted_values, unname(expected), tolerance = 1e-6)
})

test_that("Cox-Aalen combines cumulative additive and proportional components", {
    skip_if_not_installed("timereg")
    set.seed(943)
    d <- data.frame(id = 1:100, time = rexp(100), status = rbinom(100, 1, .6),
                    A = rbinom(100, 1, .5), L = rbinom(100, 1, .5))
    p <- data.frame(time = c(.5, 1), start = c(0, .2), A = c(0, 1), L = c(0, 1))
    d$status[order(d$time, decreasing = TRUE)[1:10]] <- 0
    prop <- timereg::prop
    formula <- survival::Surv(time, status) ~ -1 + factor(L) + prop(A)
    result <- learn_cox_aalen(formula, d, p, time_ref = "start")
    oracle <- timereg::cox.aalen(formula, data = d, n.sim = 0, robust = 0, resample.iid = 0)
    times <- sort(unique(c(p$start, p$time)))
    survival <- stats::predict(oracle, newdata = p, times = times, se = FALSE)$S0
    expected <- vapply(seq_len(nrow(p)), function(i) {
        survival[i, match(p$time[i], times)] / survival[i, match(p$start[i], times)]
    }, numeric(1))
    expect_s3_class(result$fit, "cox.aalen")
    expect_equal(result$predicted_values, expected, tolerance = 1e-8)

})

test_that("timereg learners reject randomized ties and invalid intervals", {
    skip_if_not_installed("timereg")
    d <- data.frame(id = 1:4, time = c(1, 1, 2, 3), status = c(1, 1, 0, 0))
    for (learner in list(learn_aalen, learn_cox_aalen)) {
        expect_error(learner("survival::Surv(time, status) ~ 1", d), "ties randomly")
        d$status <- 0
        expect_equal(learner("survival::Surv(time, status) ~ 1", d)$predicted_values, rep(1, 4))
        p <- data.frame(time = 1, start = 2)
        expect_error(learner("survival::Surv(time, status) ~ 1", d, p, time_ref = "start"), "ordered")
        expect_error(learner("survival::Surv(time, status, type = 'left') ~ 1", d), "start-stop")
        d$status <- c(1, 1, 0, 0)
    }
})

test_that("negative fitted Aalen jumps are explicitly truncated", {
    skip_if_not_installed("timereg")
    d <- data.frame(id = 1:4, time = 1:4, status = c(1, 0, 0, 0))
    fitted <- learn_aalen("survival::Surv(time, status) ~ 1", d)
    fitted$fit$cum[2, 2] <- -.2
    expect_warning(probability <- learn_aalen("ignored", d, data.frame(time = 2), reuse_fit = fitted), "truncated at zero")
    expect_equal(probability$predicted_values, 1)
})

test_that("time-varying Aalen and Cox-Aalen learners run in the public workflow", {
    skip_if_not_installed("timereg")
    set.seed(1134)
    n <- 100
    baseline <- data.frame(id = seq_len(n), A = rbinom(n, 1, .5), L = rbinom(n, 1, .5))
    last_A <- rbinom(n, 1, .5)
    measured_L <- rbinom(n, 1, .5)
    terminal <- sample(c("C", "Y", "D"), n, replace = TRUE)
    terminal[(n - 9):n] <- "Y"
    event_data <- data.frame(id = rep(seq_len(n), each = 3),
        time = as.vector(rbind(rep(.5, n), rep(1, n), 1.2 + seq_len(n) / n)),
        event = as.vector(rbind(rep("L", n), rep("A", n), terminal)),
        A = as.vector(rbind(baseline$A, last_A, last_A)),
        L = as.vector(rbind(measured_L, measured_L, rep(999, n))))
    prepare <- function(events) {
        x <- rtmle_init(c(0, 3), name_id = "id", name_time = "time", name_outcome = "Y",
                        name_competing = "D", name_censoring = "C", continuous_time = TRUE)
        x <- add_baseline_data(x, baseline)
        x <- add_continuous_data(x, event_data = events, time_covariates = "L",
                                 treatment_events = "A", covariate_events = "L")
        x <- prepare_rtmle_data(x, min_events = 0)
        x <- regime(x, "Always_A", intervention = 1, treatment_variables = "A")
        x <- target(x, "risk", regimes = "Always_A")
        model_formula(x, verbose = FALSE)
    }
    x <- prepare(event_data)
    changed <- event_data
    changed$L[changed$event %in% c("C", "Y", "D")] <- -999
    changed <- prepare(changed)
    for (learner in c("learn_aalen", "learn_cox_aalen")) {
        formula <- if (learner == "learn_aalen") ~ -1 + factor(L) else ~ -1 + factor(L) + prop(A)
        specification <- list(name = learner, fun = learner, time_varying = TRUE,
                               censoring_formula = formula)
        result <- run_rtmle(x, learner = "learn_glm", propensity_learner = "learn_glm",
                            censoring_learner = specification, estimator = "one-step",
                            save_fitted_objects = TRUE)
        expect_true(all(is.finite(result$estimate$Main_analysis$Estimate)))
        fit <- result$continuous_nuisance$Always_A$models$censoring$time_varying
        expect_false(is.null(fit$fit))
        other <- run_rtmle(changed, learner = "learn_glm", propensity_learner = "learn_glm",
                           censoring_learner = specification, estimator = "one-step",
                           save_fitted_objects = FALSE)
        expect_equal(result$estimate$Main_analysis$Estimate, other$estimate$Main_analysis$Estimate)
        expect_null(other$continuous_nuisance$Always_A$models$censoring$time_varying[["fit", exact = TRUE]])
    }
})


test_that("timereg learners reject unidentifiable additive risk sets", {
    skip_if_not_installed("timereg")
    d <- data.frame(id = c(1, 2, 2, 3, 3), start = c(0, 0, 2, 0, 2),
                    time = c(1, 2, 4, 2, 5), status = c(1, 0, 0, 0, 0), A = c(0, 0, 1, 0, 0))
    for (learner in list(learn_aalen, learn_cox_aalen)) {
        expect_error(learner("survival::Surv(start, time, status) ~ A", d, time_ref = "start"),
                     "not identifiable.*time 1")
        expect_error(learner("survival::Surv(start, time, status) ~ 1 + offset(A)", d,
                             time_ref = "start"), "offsets")
    }
})

test_that("unused factor levels do not alter additive censoring estimates", {
    skip_if_not_installed("timereg")
    d <- data.frame(id = 1:8, time = 1:8, status = c(1, 1, 1, 1, 0, 0, 0, 0),
                    group = factor(rep(c("a", "b"), 4)))
    other <- d
    other$group <- factor(other$group, levels = c("a", "b", "unused"))
    for (learner in list(learn_aalen, learn_cox_aalen)) {
        a <- learner("survival::Surv(time, status) ~ -1 + group", d)
        b <- learner("survival::Surv(time, status) ~ -1 + group", other)
        expect_equal(a$predicted_values, b$predicted_values)
    }
})

test_that("additive censoring learners reject random effects", {
    skip_if_not_installed("timereg")
    d <- data.frame(id = 1:6, time = 1:6, status = c(1, 0, 1, 0, 0, 0))
    for (learner in list(learn_aalen, learn_cox_aalen)) {
        expect_error(learner("survival::Surv(time, status) ~ frailty(id)", d), "Random effects")
        expect_error(learner("survival::Surv(time, status) ~ (1 | id)", d), "Random effects")
    }
})


test_that("timereg default fits preserve the caller's random seed", {
    skip_if_not_installed("timereg")
    set.seed(823)
    d <- data.frame(id = 1:40, time = rexp(40), status = rbinom(40, 1, .7),
                    frailty = rbinom(40, 1, .5))
    d$status[order(d$time, decreasing = TRUE)[1:8]] <- 0
    previous <- .Random.seed
    learn_aalen("survival::Surv(time, status) ~ frailty", d)
    expect_identical(.Random.seed, previous)
    learn_cox_aalen("survival::Surv(time, status) ~ prop(frailty)", d)
    expect_identical(.Random.seed, previous)
})
