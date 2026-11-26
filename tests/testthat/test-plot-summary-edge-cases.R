test_that("risk plots validate analysis tables before constructing layers", {
    data(rtmle_object)
    x <- rtmle_object
    risk_plot <- plot(x)
    expect_s3_class(risk_plot$layers[[1L]]$geom, "GeomRibbon")
    expect_error(plot(x, bootstrap_conf_int = NA), "single TRUE/FALSE")
    expect_error(plot(x, analysis = "absent"), "not found")
    expect_error(plot(x, targets = character()), "at least one target")
    expect_error(plot(x, targets = "unknown"), "Target.*not found")

    x$estimate$other <- "unfinished"
    expect_error(plot(x, analysis = "other"), "tabular")
    main <- x$estimate$Main_analysis
    x$estimate$other <- data.frame(Regime = "A")
    expect_error(plot(x, analysis = "other"), "missing required column.*Time_horizon")
    x$estimate$other <- main[0]
    expect_error(plot(x, analysis = "other"), "No estimates remain")

    x$estimate$Cheap_bootstrap <- list(Main_analysis = main)
    expect_error(plot(x, bootstrap_conf_int = TRUE), "missing required column.*Bootstrap_lower")
    bands <- data.frame(Regime = "Always_A", Time_horizon = 1,
                        Bootstrap_lower = 0.01, Bootstrap_upper = 0.09)
    x$estimate$Cheap_bootstrap <- list(Main_analysis = bands)
    expect_error(plot(x, bootstrap_conf_int = TRUE, targets = "Treatment effect"),
                 "No bootstrap confidence intervals remain")
    # Older single-target objects can omit the Target column in both tables.
    x$estimate$Main_analysis <- subset(main, Target == "Treatment effect", select = -Target)
    p <- plot(x, bootstrap_conf_int = TRUE, conf_int = FALSE)
    expect_identical(unique(as.character(p$data$Target)), "")
    expect_equal(p$layers[[1L]]$data$Bootstrap_upper, c(0, 0.09))
})

test_that("weight diagnostics validate regime selection and use the at-risk set", {
    x <- rtmle_init(0:2, name_id = "id", name_outcome = "Y")
    x$followup <- data.table::data.table(last_interval = c(0L, 1L))
    x$regimes <- list()
    expect_error(plot_IPW(x), "contains no regimes")
    x$regimes <- list(A = list())
    expect_error(plot_IPW(x), "None of the regimes has been fitted")
    # Only the second participant remains at risk at the second node.
    x$regimes$A <- list(ipw_last_nodes = 1:2,
                        intervention_last_nodes = list(node_0 = NA, node_1 = NA),
                        cumulative_intervention_probs = cbind(c(.9, .8), c(.7, .6)))
    x$prepared_data <- data.frame(id = 1:2)
    expect_error(plot_IPW(x, regimes = "B"), "Unavailable regime.*B")
    p <- plot_IPW(x, regimes = "A")
    expect_equal(p$data$used_cumprobs, c(.9, .8, .6))
})

test_that("stored nuisance coefficients retain their terms and fitted horizon", {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    x <- rtmle_init(0:1, name_id = "id", name_outcome = "Y")
    x$run_time_horizons <- 1L
    x$learner <- list(name = "ridge", fun = "learn_glmnet")
    beta <- matrix(c(-1, 0.5), dimnames = list(c("(Intercept)", "age"), "coef"))
    fit <- list(fit_summary = list(selected_beta = beta))
    fixture <- data.frame(id = 1:4, Y = c(0, 0, 0, 1))
    constant <- rtmle:::fitter(intervention_node = 0L, learner = parse_learners("learn_glmnet"),
                               formula = "Y ~ 1", data = fixture, intervened_data = fixture,
                               id_variable = "id", minority_threshold = 2L, seed = 17L,
                               save_fitted_objects = TRUE)
    x$regimes <- list(A = list())
    x$models <- list(time_0 = list(
        A = list(A_0 = fit), censoring = list(C_1 = constant),
        outcome = list(Y_1 = list(fit = list(
            A = list(sequence_time_1 = fit), B = list(sequence_time_1 = constant)
        )))
    ))
    # Unfitted models have no summary from which to extract coefficients.
    x$models$time_0$censoring$C_2 <- list(fit_summary = NULL)
    x$models$time_0$outcome$Y_1$fit$C <- list(sequence_time_1 = list(fit_summary = NULL))
    coefs <- coef(x, time_horizon = 1)
    expect_false("C_2" %in% coefs$outcome)
    expect_false("C" %in% coefs$regime)
    expect_equal(coefs[node == "A", terms], c("(Intercept)", "age"))
    expect_equal(coefs[node == "A", beta], c(-1, 0.5))
    expect_equal(coefs[node == "outcome" & terms == "age", beta], 0.5)
    expect_equal(coefs[node == "censoring" | regime == "B", .N], 0L)

    expect_error(plot_model_coefficients(x, times = 12), "No matching time")
    expect_error(plot_model_coefficients(x, time_horizon = 12), "time_horizon")
    expect_error(plot_model_coefficients(x, regime = "absent"), "regime")
    outcome_plot <- plot_model_coefficients(x, times = 0L, nodes = "A",
                                            facet_by = "time", color_by = "none",
                                            plot_style = "by_outcome")
    expect_equal(outcome_plot$data$beta, 0.5)
    expect_s3_class(outcome_plot$plot$facet, "FacetWrap")
    single_panel <- plot_model_coefficients(x, times = "time_0", nodes = "A",
                                            facet_by = "none", include_intercept = TRUE,
                                            plot_style = "by_outcome")
    expect_equal(single_panel$data$beta, c(-1, 0.5))
    expect_s3_class(single_panel$plot$facet, "FacetNull")
    # The default horizon includes the baseline node.
    invisible(capture.output(p <- plot_model_coefficients(x, manhattan_color_by = "time")))
    expect_equal(unique(p$data$time), "time_0")
    expect_true("colour" %in% names(p$mapping))
    invisible(capture.output(p <- plot_model_coefficients(x, manhattan_color_by = "none")))
    expect_false("colour" %in% names(p$mapping))
    x$models <- NULL
    expect_error(plot_model_coefficients(x), "models to be a list")
    x$learner <- list(name = "forest", fun = "learn_ranger")
    expect_error(coef(x), "Can only assess regression coefficients")
    expect_error(plot_model_coefficients(x), "Can only plot regression coefficients")
})

test_that("summaries give paired standard errors for the risk difference and ratio", {
    x <- rtmle_init(0:1, name_id = "id", name_outcome = "Y")
    expect_message(empty <- summary(x), "contains no estimates")
    expect_null(empty)
    # Risk difference = 0.4 - 0.2. The paired IC is (-0.2, 0.2).
    # For the log risk ratio, the paired IC is (-0.25, 0.25).
    x$targets <- list(risk = list(regimes = c("A", "B")))
    main <- data.table::data.table(
        Target = factor("risk"), Regime = factor(c("A", "B")),
        Target_parameter = "Risk", Estimator = "tmle", Time_horizon = 1L,
        Estimate = c(.2, .4), Lower = c(.1, .2), Upper = c(.3, .6))
    x$estimate <- list(Main_analysis = main)
    x$IC <- list(risk = list(A = list(time_horizon_1 = c(-.1, .1)),
                             B = list(time_horizon_1 = c(-.3, .3))))
    expect_error(summary(x, analysis = "unknown"), "analysis called 'unknown'")
    expect_error(summary(x, targets = "unknown"), "available_targets")
    expect_error(summary(x, reference = list(risk = "missing")), "regimes")
    s <- summary(x)
    expect_equal(s[Target_parameter == "Risk_difference", Estimate], 0.2)
    expect_equal(s[Target_parameter == "Risk_difference", Standard_error], 0.2)
    expect_equal(s[Target_parameter == "Risk_ratio", Estimate], 2)
    expect_equal(s[Target_parameter == "Risk_ratio", Standard_error], 0.25)

    # B replicates are paired by index even if the stored rows are shuffled.
    draws <- data.table::data.table(B = c(1L, 2L, 2L, 1L), Target = "risk",
                                    Regime = c("A", "A", "B", "B"), Time_horizon = 1,
                                    Bootstrap_estimate = c(.1, .2, .6, .2))
    data.table::setattr(draws, "tq", 2)
    data.table::setattr(draws, "cheap_scale", .5)
    x$estimate$Cheap_bootstrap <- list(Main_analysis = draws)
    boot <- summary(x)
    # Differences (.1, .4) and ratios (2, 3) are centered on (.2, 2).
    diff_rms <- sqrt((.1^2 + .2^2) / 2)
    log_rr_rms <- log(3 / 2) / sqrt(2)
    rd <- boot[Target_parameter == "Risk_difference"]
    rr <- boot[Target_parameter == "Risk_ratio"]
    expect_equal(c(rd$Bootstrap_lower, rd$Bootstrap_upper), .2 + c(-1, 1) * diff_rms)
    expect_equal(c(rr$Bootstrap_lower, rr$Bootstrap_upper), 2 * exp(c(-1, 1) * log_rr_rms))
    x$estimate$Cheap_bootstrap <- NULL

    # Older fits store subgroup influence curves on the subgroup estimate table.
    sub <- data.table::copy(main)
    sub[, sex := "female"]
    data.table::setattr(sub, "variable", "sex")
    subgroup_ic <- list(risk = list(A = list(time_horizon_1 = c(-.04, .04)),
                                    B = list(time_horizon_1 = c(-.14, .14))))
    data.table::setattr(sub, "IC", list(female = subgroup_ic))
    x$estimate$by_sex <- sub
    s_all <- summary(x)
    s_risk <- summary(x, targets = "risk")
    expect_equal(s_all, s_risk)
    expect_equal(s_all[sex == "female" & Target_parameter == "Risk_ratio", Estimate], 2)
    # The subgroup ICs give SEs of .10 for the RD and .15 for the log RR.
    expect_equal(s_all[sex == "female" & Target_parameter == "Risk_difference",
                       Standard_error], .10)
    expect_equal(s_all[sex == "female" & Target_parameter == "Risk_ratio",
                       Standard_error], .15)
    expect_output(print(x), "female")
    # Legacy objects could also store subgroup ICs without labels.
    sub <- data.table::copy(sub)
    sub[, sex := 1L]
    data.table::setattr(sub, "IC", stats::setNames(list(subgroup_ic), ""))
    x$estimate$by_sex <- sub
    unnamed <- summary(x, analysis = "by_sex")
    expect_equal(unnamed[sex == 1L & Target_parameter == "Risk_difference", Standard_error], .10)
    expect_equal(unnamed[sex == 1L & Target_parameter == "Risk_ratio", Standard_error], .15)
})

test_that("printing an object reports pending analysis steps and learner settings", {
    x <- rtmle_init(0:7, name_id = "id", name_outcome = "Y")
    expect_output(print(x), "0, 1, \\.\\.\\., 6")
    expect_output(print(x), "contains no regimes")
    x <- regime(x, name = "always", treatment_variables = "A", intervention = 1)
    x <- target(x, name = "risk", regimes = "always")
    expect_output(print(x), "prepare_rtmle_data")
    x$learner <- parse_learners(list(name = "my_glm", fun = "learn_glm", maxit = 25))
    expect_output(print(x), "maxit=25")
    x$learner <- parse_learners(list(folds = 3, ensemble_method = "discrete",
                                     learners = list("learn_glm", "learn_glmnet")))
    expect_output(print(x), "SuperLearner \\(3-fold CV\\)")
})

test_that("a one-horizon stratified contrast uses pointwise confidence intervals", {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    data(rtmle_object)
    x <- target(rtmle_object, name = "A only", regimes = "Always_A")
    fits <- stratify(x, sex, time_horizon = 1L,
                     targets = c("A only", "Treatment effect"))
    expect_error(plot(fits, targets = "A only", target_parameter = "Risk_difference"),
                 "No estimates are available")
    for (stratum in names(fits)) {
        fits[[stratum]]$time_grid_labels <- NULL
        fits[[stratum]]$time_unit <- "Months"
    }
    p <- plot(fits, targets = "Treatment effect", target_parameter = "Risk_difference")
    expect_identical(p$labels$x, "Months")
    expect_identical(p$scales$get_scales("x")$labels, 1L)
    expect_true(any(vapply(p$layers, function(layer) inherits(layer$geom, "GeomLinerange"),
                           logical(1))))
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("plots explain which drawing dependency is missing", {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    data(rtmle_object)
    risk_plot <- without_package(rtmle:::autoplot.rtmle, "ggplot2")
    trajectories <- without_package(plot_treatment_trajectories, "ggplot2")
    coefficients <- without_package(plot_model_coefficients, "ggplot2")
    expect_error(risk_plot(rtmle_object), "Package 'ggplot2' must be installed")
    expect_error(trajectories(rtmle_object), "Package 'ggplot2' must be installed")
    expect_error(coefficients(rtmle_object), "Package 'ggplot2' is required")

    coefficients <- without_package(plot_model_coefficients, "Matrix")
    expect_error(coefficients(rtmle_object), "Package 'Matrix' is required")
    coefficients <- without_package(plot_model_coefficients, "plotly")
    invisible(capture.output(p <- coefficients(rtmle_object, times = 0)))
    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$time == "time_0"))
})
