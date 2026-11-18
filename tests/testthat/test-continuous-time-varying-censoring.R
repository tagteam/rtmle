time_varying_censoring_fixture <- function() {
    set.seed(214)
    n <- 80L
    baseline <- data.frame(id = seq_len(n), A_0 = rbinom(n, 1, .5),
                           L_0 = rnorm(n), age = rnorm(n, 50, 8))
    a1 <- rbinom(n, 1, .5)
    a3 <- rbinom(n, 1, .5)
    l1 <- rnorm(n)
    l2 <- rnorm(n)
    terminal <- rep(c("C", "Y", "D", "tauend"), length.out = n)
    last <- 4 + seq_len(n) / 40
    last[terminal == "tauend"] <- 8
    events <- rbind(data.frame(id = baseline$id, time = 1, event = "A", A = a1, L = l1),
                    data.frame(id = baseline$id, time = 2, event = "L", A = a1, L = l2),
                    data.frame(id = baseline$id, time = 3, event = "A", A = a3, L = l2),
                    data.frame(id = baseline$id, time = last, event = terminal, A = a3, L = 999))
    # Construct the reference counting-process data independently of preparation.
    oracle <- rbind(data.frame(id = baseline$id, start = 0, stop = 1, censored = 0,
                              A = baseline$A_0, L = baseline$L_0, age = baseline$age),
                    data.frame(id = baseline$id, start = 1, stop = 2, censored = 0,
                              A = a1, L = l1, age = baseline$age),
                    data.frame(id = baseline$id, start = 2, stop = 3, censored = 0,
                              A = a1, L = l2, age = baseline$age),
                    data.frame(id = baseline$id, start = 3, stop = last,
                              censored = as.integer(terminal == "C"), A = a3, L = l2, age = baseline$age))
    list(data = list(baseline_data = baseline, timevarying_data = events), oracle = oracle)
}

test_that("time-varying Cox uses start-of-interval covariates and a single fit", {
    fixture <- time_varying_censoring_fixture()
    p <- prepare_continuous_data(fixture$data, c(1.5, 7), c("A", "L"),
                                 c("A_0", "L_0", "age"), min_events = 0)
    learner <- list(name = "coxph", fun = "learn_coxph", time_varying = TRUE,
                    censoring_formula = ~ A + L + age)
    actual <- continuous_propensity(p, censoring_learner = learner)
    oracle <- survival::coxph(survival::Surv(start, stop, censored) ~ A + L + age,
                              data = fixture$oracle, x = TRUE)
    expect_equal(unname(coef(actual$models$censoring$time_varying$fit)), unname(coef(oracle)), tolerance = 1e-10)
    hazard <- survival::basehaz(oracle, centered = FALSE)
    # Directly sum baseline hazard jumps with the predictor recorded before each jump.
    expected <- vapply(p$wide_data$id, function(id) {
        rows <- fixture$oracle[fixture$oracle$id == id, ]
        endpoint <- max(rows$stop)
        jumps <- hazard$time[hazard$time < endpoint]
        increments <- diff(c(0, hazard$hazard))[hazard$time < endpoint]
        risk <- vapply(jumps, function(u) {
            row <- rows[rows$start < u & rows$stop >= u, ]
            exp(as.numeric(stats::predict(oracle, newdata = row, type = "lp", reference = "zero")))
        }, numeric(1))
        exp(-sum(risk * increments))
    }, numeric(1))
    expect_equal(actual$data$survival_censoring_pooled_1, expected, tolerance = 1e-10)

    changed <- fixture$data
    terminal <- changed$timevarying_data$event %in% c("C", "Y", "D", "tauend")
    changed$timevarying_data$L[terminal] <- -12345
    p2 <- prepare_continuous_data(changed, c(1.5, 7), c("A", "L"),
                                  c("A_0", "L_0", "age"), min_events = 0)
    expect_equal(continuous_propensity(p2, censoring_learner = learner)$data$survival_censoring_pooled_1,
                 actual$data$survival_censoring_pooled_1)
    # Pooling away later events must not remove their records from hazard fitting.
    pooled <- prepare_continuous_data(fixture$data, 7, c("A", "L"),
                                      c("A_0", "L_0", "age"), last_non_terminal_event = 0)
    expect_equal(continuous_propensity(pooled, censoring_learner = learner)$data$survival_censoring_pooled_1,
                 expected, tolerance = 1e-10)
})

test_that("pooled survival includes hazard jumps at interior recorded boundaries", {
    fixture <- time_varying_censoring_fixture()
    p <- prepare_continuous_data(fixture$data, 7, c("A", "L"), c("A_0", "L_0", "age"), min_events = 0)
    fake <- function(character_formula, data, intervened_data, time_variable, time_ref,
                      left_limit = TRUE, ...) {
        a <- intervened_data[[time_ref]]
        b <- intervened_data[[time_variable]]
        # A known hazard jump at time 2, tied with all subjects' L events.
        counted <- a < 2 & if (left_limit) b > 2 else b >= 2
        list(predicted_values = exp(-as.numeric(counted)))
    }
    specification <- continuous_censoring_specification(list(name = "known", fun = fake, time_varying = TRUE))
    call <- function(args) do.call(fake, args)
    model <- continuous_fit_censoring(p, specification, NULL, call)
    rows <- data.frame(id = 1:2, begin = c(0, 0), end = c(2, 4))
    expect_equal(continuous_predict_censoring(model, rows, "end", "begin", specification, call), c(1, exp(-1)))
    rows$end[2] <- 20
    expect_error(continuous_predict_censoring(model, rows, "end", "begin", specification, call), "beyond")
})

test_that("time-varying censoring specifications reject incompatible options", {
    spec <- function(...) continuous_censoring_specification(list(name = "cox", fun = "learn_coxph", ...))
    expect_error(spec(time_varying = NA), "time_varying")
    expect_error(spec(time_varying = TRUE, marginal_censoring = TRUE), "cannot be combined")
    expect_error(spec(censoring_formula = ~ L), "requires time_varying")
    expect_error(spec(time_varying = TRUE, censoring_formula = "L"), "one-sided")
    fixture <- time_varying_censoring_fixture()
    p <- prepare_continuous_data(fixture$data, 7, c("A", "L"), c("A_0", "L_0", "age"), min_events = 0)
    expect_error(continuous_propensity(p, censoring_learner = list(name="cox", fun="learn_coxph",
        time_varying=TRUE, censoring_formula=~L_10)), "baseline or current")
    expect_error(continuous_propensity(p, censoring_learner = list(name="cox", fun="learn_coxph",
        time_varying=TRUE), model_options=list(exclusion_rules=list("C_1$"="L"))), "same predictors")
})

test_that("the shared hazard is fitted once and factor initial values survive", {
    fixture <- time_varying_censoring_fixture()
    p <- prepare_continuous_data(fixture$data, 7, c("A", "L"), c("A_0", "L_0", "age"), min_events = 0)
    calls <- 0L
    spy <- function(..., reuse_fit = NULL) {
        if (is.null(reuse_fit)) calls <<- calls + 1L
        learn_coxph(..., reuse_fit = reuse_fit)
    }
    continuous_propensity(p, censoring_learner=list(name="spy", fun=spy,
        time_varying=TRUE, censoring_formula=~A+L))
    expect_equal(calls, 1L)
    p$baseline_data$L_0 <- factor(rep("initial", nrow(p$baseline_data)))
    p$event_data$L <- factor(rep("later", nrow(p$event_data)))
    intervals <- continuous_censoring_intervals(p)
    expect_false(anyNA(intervals$L))
    expect_true(all(as.character(intervals$L[intervals$.rtmle_start == 0]) == "initial"))
})

test_that("the public Cox workflow constructs interval data from ordinary event inputs", {
    fixture <- time_varying_censoring_fixture()
    x <- continuous_analysis(fixture$data, horizons = 7, min_events = 0)
    result <- run_rtmle(x, estimator = "one-step", save_fitted_objects = TRUE,
        censoring_learner = list(name = "cox", fun = "learn_coxph", time_varying = TRUE,
                                 censoring_formula = ~ A + L + age))
    expect_true(all(is.finite(result$estimate$Main_analysis$Estimate)))
    fitted <- result$continuous_nuisance$Always_A$models$censoring$time_varying$fit
    expect_identical(attr(fitted$y, "type"), "counting")
    expect_equal(nrow(fitted$y), nrow(fixture$data$timevarying_data))
})
