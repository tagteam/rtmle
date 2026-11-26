make_one_interval_trial <- function() {
    id <- seq_len(80)
    treatment <- rep(0:1, each = 40)
    # Ten events among the 40 untreated, 30 among the 40 treated.
    outcome <- c(rep(0:1, c(30, 10)), rep(0:1, c(10, 30)))
    x <- rtmle_init(0:1, name_id = "id", name_outcome = "Y")
    x <- add_wide_data(x, outcome_data = data.frame(id, Y_1 = outcome),
                       timevar_data = list(A = data.frame(id, A_0 = treatment)))
    x <- prepare_rtmle_data(x)
    x <- regime(x, "untreated", treatment_variables = "A", intervention = 0)
    x <- regime(x, "treated", treatment_variables = "A", intervention = 1)
    x <- target(x, "risk", c("untreated", "treated"))
    model_formula(x, verbose = FALSE)
}

test_that("one-interval TMLE recovers the arm risks and their influence curves", {
    x <- make_one_interval_trial()
    utils::capture.output(
        messages <- utils::capture.output(
            fit <- run_rtmle(x, save_fitted_objects = TRUE,
                             verbose = TRUE, progressbar = 1),
            type = "message"
        )
    )
    expect_true(any(grepl("Running target: risk", messages, fixed = TRUE)))
    expect_true(any(grepl("Fitting sequential regression", messages, fixed = TRUE)))
    estimates <- fit$estimate$Main_analysis
    outcome <- x$prepared_data$Y_1
    for (arm in c("untreated", "treated")) {
        treated <- arm == "treated"
        risk <- if (treated) 0.75 else 0.25
        # Balanced allocation gives P(A = a) = 1/2.
        ic <- 2 * (as.integer(x$prepared_data$A_0 == "1") == treated) * (outcome - risk)
        row <- estimates[as.character(estimates$Regime) == arm, ]
        expect_equal(row$Estimate, risk, tolerance = 1e-8)
        expect_equal(row$Standard_error, stats::sd(ic) / sqrt(80), tolerance = 1e-8)
        expect_equal(fit$IC$risk[[arm]]$time_horizon_1, ic, tolerance = 1e-8)
        outcome_fit <- fit$models$time_0$outcome$Y_1$fit[[arm]]$sequence_time_1$fit
        newdata <- data.frame(A_0 = factor(0:1))
        expect_equal(as.numeric(stats::predict(outcome_fit, newdata, type = "response")),
                     c(.25, .75),
                     tolerance = 1e-7)
    }

    gformula <- run_rtmle(x, estimator = "g-formula")
    estimates <- gformula$estimate$Main_analysis
    expect_equal(estimates$Estimate[match(c("untreated", "treated"), estimates$Regime)],
                 c(0.25, 0.75), tolerance = 1e-8)
})

test_that("initialization rejects invalid grids, axis labels and censoring levels", {
    expect_error(rtmle_init(c(1, 2), "id", name_outcome = "Y"), "starting with 0")
    expect_error(rtmle_init(0, "id", name_outcome = "Y"), "length > 1")
    expect_error(rtmle_init(0:1, "id", name_outcome = "Y", time_grid_labels = "start"),
                 "same length")
    expect_error(rtmle_init(0:1, "id", name_outcome = "Y", name_censoring = "C",
                            censored_levels = c("lost", "observed"), censored_label = "lost"),
                 "second value")
})

test_that("fitting requires models, valid targets and matching prepared data", {
    x <- rtmle_init(0:1, "id", name_outcome = "Y")
    expect_error(run_rtmle(x), "model_formula")
    expect_error(target(x, "risk", "absent"), "regimes are not defined")
    x <- make_one_interval_trial()
    expect_error(run_rtmle(x, targets = "unknown"), "Requested targets")
    expect_error(run_rtmle(x, time_horizon = 2), "time_horizon")
    expect_error(run_rtmle(x, time_horizon = 0), "time_horizon")
    expect_error(run_rtmle(x, keep_influence = TRUE), "Use subset")
    x$targets$incomplete <- list(regimes = "absent")
    expect_error(run_rtmle(x), "Referenced by target 'incomplete': absent")
    x$regimes <- NULL
    expect_error(run_rtmle(x), "Available regimes: [(]none[)]")
    x$targets <- NULL
    expect_error(run_rtmle(x), "contains no targets")

    x <- make_one_interval_trial()
    x$prepared_data$id <- NULL
    expect_error(run_rtmle(x), "Cannot see id variable id")
})

test_that("the fluctuation update solves the weighted intercept equation", {
    update <- rtmle:::tmle_update
    y <- c(0, 0, 1, 1, 1, NA_real_)
    probability <- c(0.5, 1, 0.25, 1, 1, NA_real_)
    at_risk <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
    matched <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
    # The three included outcomes have inverse probability weights 2, 1 and 4.
    fitted <- update(y, rep(0, 6), probability, at_risk, matched, k = 1, regime = "A")
    expect_equal(as.numeric(fitted), rep(4 / 7, 6), tolerance = 1e-8)
    unweighted <- update(y, rep(0, 6), NULL, at_risk, matched, k = 1, regime = "A")
    expect_equal(as.numeric(unweighted), rep(1 / 3, 6), tolerance = 1e-8)

    offsets <- stats::qlogis(seq(0.1, 0.6, by = 0.1))
    no_match <- update(y, offsets, probability, at_risk, rep(FALSE, 6), k = 2, regime = "B")
    expect_equal(as.numeric(no_match), stats::plogis(offsets))
    diagnostic <- attr(no_match, "diagnostics")
    expect_equal(diagnostic$Step, 2)
    expect_equal(diagnostic$Regime, "B")
    expect_equal(diagnostic$Event, "No positive weights")
    probability[1] <- NA_real_
    expect_error(update(y, offsets, probability, at_risk, matched, k = 1, regime = "A"),
                 "NA in weights")
})

test_that("clearing a partial fit drops its pseudo-outcome and retains original outcomes", {
    x <- make_one_interval_trial()
    original <- data.table::copy(x$prepared_data)
    x$prepared_data$rtmle_predicted_outcome <- rep(0.5, nrow(original))
    clean <- rtmle:::clear_rtmle_fit(x)
    expect_equal(clean$prepared_data, original)
    expect_equal(clean$models$time_0$outcome$Y_1$formula,
                 x$models$time_0$outcome$Y_1$formula)
})

test_that("estimation rejects zero probabilities and unusable fluctuation offsets", {
    x <- make_one_interval_trial()
    expect_error(run_rtmle(x, weight_truncation = c(0, 0)),
                 "Exactly zero intervention probabilities")
    # With no events and clipping disabled, a zero prediction gives an infinite offset.
    x$prepared_data$Y_1[] <- 0
    utils::capture.output(
        expect_error(run_rtmle(x, prediction_range = c(0, 1)),
                     "Fluctuation model.*time node: 1"),
        type = "message"
    )
})

test_that("a regime can start after baseline when censoring is modeled", {
    id <- seq_len(80)
    # Within every observed treatment history exactly one quarter have an event.
    treatments <- data.frame(id, A_0 = rep(rep(0:1, each = 4), 10),
                             A_1 = rep(0:1, each = 40))
    outcomes <- data.frame(id, Y_1 = 0, Y_2 = rep(c(0, 0, 0, 1), 20),
                           C_1 = "uncensored", C_2 = "uncensored")
    x <- rtmle_init(0:2, "id", name_outcome = "Y", name_censoring = "C")
    x <- add_wide_data(x, outcome_data = outcomes, timevar_data = list(A = treatments))
    x <- prepare_rtmle_data(x)
    x <- regime(x, "start_at_1", intervention = data.frame(time_node = 1,
                                                           A = factor(1, levels = 0:1)),
                expand = FALSE)
    x <- target(x, "risk", "start_at_1")
    x <- model_formula(x, verbose = FALSE)
    fit <- run_rtmle(x)
    expect_equal(fit$estimate$Main_analysis$Estimate, 0.25, tolerance = 1e-8)
    expect_true(is.na(fit$regimes$start_at_1$intervention_last_nodes[["node_0"]]))
})

test_that("fitting records each regime with no adherent participants at the second visit", {
    id <- seq_len(80)
    a0 <- rep(rep(0:1, each = 4), 10)
    # Everyone switches treatment at the second visit.
    x <- rtmle_init(0:2, "id", name_outcome = "Y")
    x <- add_wide_data(x,
        outcome_data = data.frame(id, Y_1 = 0, Y_2 = rep(c(0, 0, 0, 1), 20)),
        timevar_data = list(A = data.frame(id, A_0 = a0, A_1 = 1 - a0)))
    x <- add_baseline_data(x, data.frame(id, Z = rep(c(-1, 1), each = 40)))
    x <- prepare_rtmle_data(x)
    x <- regime(x, "always", treatment_variables = "A", intervention = 1)
    x <- regime(x, "never", treatment_variables = "A", intervention = 0)
    x <- target(x, "risk", c("always", "never"))
    # A_0 and A_1 are perfectly collinear in these data.
    x <- model_formula(x, exclusion_rules = list(A_1 = "A_0", Y_2 = "A_0"), verbose = FALSE)
    fit <- run_rtmle(x)
    diagnostic <- fit$diagnostics$no_positive_weights
    expect_setequal(diagnostic$Regime, c("always", "never"))
    expect_equal(diagnostic$Step, c(2, 2))
    expect_equal(diagnostic$Event, rep("No positive weights", 2))
})
