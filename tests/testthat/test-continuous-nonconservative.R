test_that("the non-conservative update adds each subject's censoring contribution", {
    treatment <- c(1, 1, 1, 0, 0)
    event <- c("C", "D", "Y", "Y", "tauend")
    d <- list(baseline_data = data.frame(id = 1:5, A_0 = treatment, L_0 = 0),
              timevarying_data = data.frame(id = 1:5, time = c(1, 1, 2, 1, 3),
                                            event = event, A = treatment, L = 0))
    x <- continuous_analysis(d, horizons = c(1, 2), regimes = "Always_A")
    fit_mean <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = rep(mean(data$rtmle_outcome), nrow(intervened_data)))
    }
    fit <- function(x, conservative, time_horizon, grid = 4L) {
        run_rtmle(x, estimator = "one-step", learner = fit_mean, propensity_learner = "glm",
                  censoring_learner = list(name = "cox", fun = "learn_coxph",
                                           marginal_censoring_formula = ~ 1),
                  conservative = conservative, integration_grid_size = grid,
                  time_horizon = time_horizon, refit = FALSE)
    }
    initial <- fit(x, TRUE, 1:2)
    updated <- fit(initial, FALSE, 2)

    # At time 1 there is one dropout, one outcome, and one competing event.
    # The intercept-only Cox model gives G(1) = exp(-1/5). The integrand
    # uses q(1) = 1/5 after the outcome jump and S(1-) = 1 before any event.
    G <- exp(-1/5)
    q_at_jump <- 1/5
    q_tau <- (1 + 1/G)/5
    weight <- treatment / mean(treatment)
    phi_c <- weight * (q_tau - q_at_jump) / G * ((event == "C") - (1 - G))
    phi_d <- initial$IC$risk$Always_A$time_horizon_2
    e <- updated$estimate$Main_analysis
    expect_equal(e$Time_horizon, 2L)
    expect_equal(e$ICE_IPCW, q_tau)
    expect_equal(updated$IC$risk$Always_A$time_horizon_2, phi_d + phi_c)
    expect_equal(e$Estimate, q_tau + mean(phi_d + phi_c))
    expect_equal(e$Standard_error, sd(phi_d + phi_c)/sqrt(5))
    # With five cells the three simultaneous jumps lie inside a cell.
    off_grid <- fit(initial, FALSE, 2, grid = 5L)
    expect_equal(off_grid$IC$risk$Always_A$time_horizon_2, phi_d + phi_c)

    # Keep results from matching fits; changing the quadrature invalidates
    # estimates from the earlier grid even when refit = FALSE.
    both <- fit(updated, FALSE, 1)
    expect_equal(both$estimate$Main_analysis$Time_horizon, 1:2)
    changed <- fit(both, FALSE, 2, grid = 8L)
    expect_equal(changed$estimate$Main_analysis$Time_horizon, 2L)
    expect_equal(changed$IC, updated$IC)
})

test_that("administrative follow-up without dropout gives zero censoring correction", {
    x <- continuous_analysis()
    original <- run_rtmle(x, estimator = "one-step")
    updated <- run_rtmle(x, estimator = "one-step", conservative = FALSE,
                         integration_grid_size = 4L)
    expect_equal(updated$estimate, original$estimate)
    expect_equal(updated$IC, original$IC)
    expect_error(run_rtmle(x, estimator = "one-step", conservative = FALSE,
                           integration_grid_size = 1), "integration_grid_size")
    expect_error(run_rtmle(x, estimator = "one-step", conservative = FALSE,
                           event_learner = list(name = "cox", fun = "learn_coxph",
                                                marginal_censoring = TRUE)),
                 "marginal|conditional")
})

test_that("quadrature uses the pooled endpoint, interval history and adherence weight", {
    # The first subject drops out inside a cell. The second deviates from the
    # regime at its first event. The third drops out after the horizon.
    grid <- 0:4
    entry <- c(1.25, 2, 0.5)
    data <- data.table::data.table(
        id = 1:3, event_0 = "A", time_0 = 0, A_0 = 1, propensity_0 = 0.8,
        event_1 = c("L", "A", "L"), time_1 = entry, A_1 = c(1, 0, 1),
        propensity_1 = c(NA, 0.7, NA), survival_censoring_1 = exp(-0.2 * entry),
        event_2 = "A", time_2 = c(2, 2.1, 1),
        event_pooled_2 = "C", time_pooled_2 = c(2.5, 2.25, 5))
    nuisance <- list(data = data, prepared_data_object = list(
        info = data.frame(time_horizon = 4, last_event = 2L)))
    q <- list(matrix(0, nrow = 3, ncol = length(grid)),
              t(vapply(entry, function(t) pmax(grid - t, 0) / 4, numeric(length(grid)))))
    local_mocked_bindings(
        continuous_grid_survival = function(nuisance, k, pooled, rows, u, from, jump_rows) {
            start <- nuisance$data[[paste0("time_", k - 1L)]][rows]
            list(event = exp(-0.1 * (from - start)),
                 event_at_jump = exp(-0.1 * (u - start)[jump_rows]),
                 censoring = exp(-0.2 * (u - start)))
        }, .package = "rtmle")

    weights <- continuous_weights(data, 1, 2)
    result <- continuous_censoring_correction(nuisance, 4, grid, q, weights)
    contribution <- function(start, right, dropout) {
        width <- diff(c(start, right))
        q_difference <- (4 - ceiling(right)) / 4
        G <- exp(-0.2 * (right - start))
        S_right <- exp(-0.1 * (right - start))
        S_left <- exp(-0.1 * (right - width - start))
        jump <- (right == dropout) / S_right
        compensator <- (1 - exp(-0.2 * width)) / S_left
        weight <- 1 / (0.8 * exp(-0.2 * start))
        weight * sum(q_difference / G * (jump - compensator))
    }
    expect_equal(result, c(contribution(1.25, c(2, 2.5), 2.5), 0,
                           contribution(0.5, 1:4, 5)))
})

test_that("event survival uses the next event and each interval's preceding history", {
    baseline <- data.frame(id = 1:8, A_0 = rep(0:1, 4), L_0 = 0)
    events <- data.frame(
        id = c(1:8, 1:4),
        time = c(0.8, 0.9, 1, 1.1, 1, 1, 1, 3, 2, 2, 2, 3),
        event = c("A", "L", "A", "L", "C", "Y", "D", "tauend", "Y", "C", "D", "tauend"),
        A = c(1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1),
        L = c(0, 3, 0, 1, 0, 0, 0, 0, 0, 3, 0, 1))
    prepared <- prepare_continuous_data(list(baseline_data = baseline, timevarying_data = events),
        c(1, 2), c("A", "L"), c("A_0", "L_0"), min_events = 0)

    calls <- list()
    record <- function(character_formula, data, intervened_data, time_variable, ...) {
        formula <- stats::as.formula(character_formula)
        response <- stats::model.response(stats::model.frame(formula, data))
        calls[[time_variable]] <<- list(
            ids = data$id, duration = response[, "time"], status = response[, "status"],
            predictors = all.vars(formula[[3L]]))
        list(predicted_values = rep(1, nrow(intervened_data)))
    }
    constant <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = rep(0.5, nrow(intervened_data)))
    }
    continuous_propensity(prepared, learner = constant,
                          censoring_learner = list(fun = constant, name = "constant",
                                                   marginal_censoring = FALSE),
                          event_learner = record)
    # The ordinary first interval ends at A, L, Y or D. When later updates
    # are pooled into the terminal record, only Y and D are events.
    expect_equal(unname(calls$time_1$status), c(1, 1, 1, 1, 0, 1, 1, 0))
    expect_equal(unname(calls$time_pooled_1$status), c(1, 0, 1, 0, 0, 1, 1, 0))
    second <- calls$time_pooled_2
    expect_equal(second$ids, 1:4)
    expect_equal(unname(second$duration), c(1.2, 1.1, 1, 1.9))
    expect_true(all(c("time_1", "A_1", "L_1") %in% second$predictors))
    expect_false(any(c("A_2", "L_2") %in% second$predictors))
})
