test_that("grid regressions truncate the next event and retain its continuation at tau", {
    # Subjects 1 and 2 have a visit at 0.7. Only subject 1 later has an
    # outcome. Subject 3 has an outcome without a visit. Tau is 2.
    d <- list(
        baseline_data = data.frame(id = 1:4, A_0 = 1),
        timevarying_data = data.frame(
            id = c(1, 1, 2, 2, 3, 4), time = c(0.7, 1.7, 0.7, 2, 1.5, 2),
            event = c("L", "Y", "L", "tauend", "Y", "tauend"), A = 1
        )
    )
    prepared <- prepare_continuous_data(d, 2, "A", "A_0", min_events = 0)
    nuisance <- continuous_propensity(prepared)
    fit <- function(grid = NULL) {
        continuous_sequential_regression(
            nuisance, th = 2, a = 1, learner = continuous_nuisance_learner("learn_glm"),
            estimator = "one-step", lag = NULL, seed = 1, prediction_range = c(0.0001, 0.9999),
            stratify = FALSE, intervene_all = FALSE,
            intervention_hook = NULL, save_fitted_objects = FALSE, q_grid = grid
        )
    }
    grid_fit <- fit(c(0, 0.5, 0.7, 1, 2))
    horizon_fit <- fit()
    expect_identical(grid_fit[names(horizon_fit)], horizon_fit)
    # q0(1) = (2/4) * (1/2), even though no outcome occurs by time 1.
    expect_equal(unname(grid_fit$q_grid[[1]]),
                 matrix(rep(c(0, 0, 0.25, 0.25, 0.5), each = 4), nrow = 4))
    q1 <- grid_fit$q_grid[[2]]
    expect_equal(q1[, 1:4], matrix(0, nrow = 4, ncol = 4))
    expect_equal(q1[, 5], c(0.5, 0.5, 0, 0))
})

test_that("auxiliary grid fits preserve the horizon fits of a stochastic learner", {
    d <- continuous_example()
    events <- d$timevarying_data
    events$time[events$event == "Y"] <- 1.4
    visits <- transform(events, time = 0.4, event = "L")
    d$timevarying_data <- rbind(visits, events)
    prepared <- prepare_continuous_data(d, 2, c("A", "L"), c("A_0", "L_0"), min_events = 0)
    nuisance <- continuous_propensity(prepared)
    random_learner <- function(character_formula, data, intervened_data, ...) {
        list(predicted_values = stats::runif(nrow(intervened_data)))
    }
    fit <- function(grid, grid_learner = random_learner) {
        set.seed(101)
        result <- continuous_sequential_regression(
            nuisance, th = 2, a = 1, learner = continuous_nuisance_learner(random_learner),
            estimator = "one-step", lag = NULL, seed = NULL, prediction_range = c(0.0001, 0.9999),
            stratify = FALSE, intervene_all = FALSE,
            intervention_hook = NULL, save_fitted_objects = FALSE, q_grid = grid,
            grid_learner = continuous_nuisance_learner(grid_learner)
        )
        list(estimate = result$estimate, ice_ipcw = result$ice_ipcw, ic = result$ic,
             next_draw = stats::runif(1))
    }
    expected <- fit(NULL)
    expect_identical(fit(seq(0, 2, by = 0.5)), expected)
    expect_identical(fit(seq(0, 2, by = 0.25)), expected)
    invalid_grid_learner <- function(...) list(predicted_values = NA_real_)
    expect_error(fit(seq(0, 2, by = 0.5), invalid_grid_learner),
                 "Grid outcome regression failed at event 2, u = 1.5 .*Outcome learner")
})
