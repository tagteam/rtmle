counting_cox_fixture <- function() {
    n <- 24L
    data.frame(id = rep(seq_len(n), each = 2),
        start = rep(c(0, 2), n), stop = as.vector(rbind(rep(2, n), 3 + seq_len(n) / 3)),
        status = as.vector(rbind(rep(0, n), as.integer(seq_len(n) %% 3 != 0))),
        Z = sin(seq_len(2 * n)), group = rep(rep(c("a", "b"), n / 2), each = 2))
}

test_that("counting-process Cox survival agrees with a separate absolute-time fit", {
    d <- counting_cox_fixture()
    p <- data.frame(start = c(0, 2, 4), stop = c(5, 7, 8), Z = c(-.5, 0, .5))
    direct <- survival::coxph(survival::Surv(start, stop, status) ~ Z, data = d)
    hazard <- survival::basehaz(direct, centered = FALSE)
    increments <- diff(c(0, hazard$hazard))
    expected <- vapply(seq_len(nrow(p)), function(i) {
        selected <- hazard$time > p$start[i] & hazard$time < p$stop[i]
        exp(-sum(increments[selected]) * exp(stats::coef(direct)[[1]] * p$Z[i]))
    }, numeric(1))
    actual <- learn_coxph("survival::Surv(start, stop, status) ~ Z", d, p,
        time_variable = "stop", time_ref = "start")
    expect_equal(stats::coef(actual$fit), stats::coef(direct), tolerance = 1e-12)
    expect_equal(actual$predicted_values, expected, tolerance = 1e-12)

    reused <- learn_coxph("unused", d, p, time_variable = "stop", time_ref = "start",
                          reuse_fit = actual)
    expect_equal(reused$predicted_values, expected, tolerance = 1e-12)
})

test_that("counting-process Cox handles strata and internal interval boundaries", {
    d <- counting_cox_fixture()
    p <- data.frame(start = c(0, 2), stop = c(5, 7), Z = c(-.5, .5), group = c("a", "b"))
    direct <- survival::coxph(survival::Surv(start, stop, status) ~ Z + strata(group), data = d)
    hazard <- survival::basehaz(direct, centered = FALSE)
    expected <- vapply(seq_len(nrow(p)), function(i) {
        h <- hazard[as.character(hazard$strata) == p$group[i], ]
        jumps <- diff(c(0, h$hazard))
        exp(-sum(jumps[h$time > p$start[i] & h$time < p$stop[i]]) *
            exp(stats::coef(direct)[[1]] * p$Z[i]))
    }, numeric(1))
    actual <- learn_coxph("survival::Surv(start, stop, status) ~ Z + strata(group)", d, p,
        time_variable = "stop", time_ref = "start")
    expect_equal(actual$predicted_values, expected, tolerance = 1e-12)

    # A censoring jump at a shared interior boundary belongs to the interval
    # ending there. The final endpoint alone uses a left limit.
    d <- data.frame(start = 0, stop = c(2, 3, 4), status = c(1, 1, 0))
    p <- data.frame(start = 0, stop = 2)
    left <- learn_coxph("survival::Surv(start, stop, status) ~ 1", d, p,
        time_variable = "stop", time_ref = "start")
    right <- learn_coxph("unused", d, p, time_variable = "stop", time_ref = "start",
        reuse_fit = left, left_limit = FALSE)
    expect_equal(left$predicted_values, 1)
    expect_equal(right$predicted_values, exp(-1 / 3), tolerance = 1e-12)
})

test_that("counting-process Cox checks times and detects censoring in the status column", {
    d <- counting_cox_fixture()
    d$status <- 0
    # The start column contains 1 here, but this is not a censoring event.
    d$start[d$start == 2] <- 1
    result <- learn_coxph("survival::Surv(start, stop, status) ~ Z", d,
        time_variable = "stop", time_ref = "start")
    expect_null(result$fit)
    expect_equal(result$predicted_values, rep(1, nrow(d)))
    bad <- d
    bad$start[1] <- -1
    expect_error(learn_coxph("survival::Surv(start, stop, status) ~ Z", bad,
        time_variable = "stop", time_ref = "start"), "nonnegative")
    bad <- d
    bad$stop[1] <- Inf
    expect_error(learn_coxph("survival::Surv(start, stop, status) ~ Z", bad,
        time_variable = "stop", time_ref = "start"), "finite")
    bad <- d
    bad$stop[1] <- bad$start[1]
    expect_error(suppressWarnings(learn_coxph(
        "survival::Surv(start, stop, status) ~ Z", bad,
        time_variable = "stop", time_ref = "start")), "missing|NA")
    bad <- d
    bad$status[1] <- NA_integer_
    expect_error(learn_coxph("survival::Surv(start, stop, status) ~ Z", bad,
        time_variable = "stop", time_ref = "start"), "missing|NA")
    bad <- d
    bad$Z[1] <- NA_real_
    expect_error(learn_coxph("survival::Surv(start, stop, status) ~ Z", bad,
        time_variable = "stop", time_ref = "start"), "missing|NA")
    expect_error(learn_coxph("unused", d, reuse_fit = result, left_limit = NA), "left_limit")
})

test_that("Cox censoring rejects random effects while retaining repeated interval IDs", {
    d <- counting_cox_fixture()
    for (term in c("frailty(id)", "frailty.gamma(id)", "survival::frailty(id)",
                   "survival::frailty.gaussian(id)", "(1 | id)", "(1 || id)")) {
        expect_error(learn_coxph(paste("survival::Surv(start, stop, status) ~ Z +", term),
            d, time_variable = "stop", time_ref = "start"), "frailty or random effects")
    }
    d$frailty <- d$Z
    ordinary <- learn_coxph("survival::Surv(start, stop, status) ~ frailty", d,
        time_variable = "stop", time_ref = "start")
    expect_true(all(is.finite(ordinary$predicted_values)))
    expect_equal(nrow(ordinary$fit$y), nrow(d))
})
