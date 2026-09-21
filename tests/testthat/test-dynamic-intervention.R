library(testthat)
library(rtmle)
library(data.table)

make_dynamic_intervention_fixture <- function(repeats = 1L) {
    histories <- data.table(
        L_0 = c(1, 1, 0, 0, 0, 0, 0, 0),
        A_0 = c(0, 1, 1, 1, 1, 0, 1, 0),
        L_1 = c(1, 1, 1, 1, 0, 0, 0, 0),
        A_1 = c(0, 0, 0, 1, 1, 1, 0, 0)
    )
    d <- histories[rep(seq_len(NROW(histories)), times = repeats)]
    n <- NROW(d)
    d[, `:=`(
        id = seq_len(n),
        W = seq(-1, 1, length.out = n),
        Y_1 = as.integer(seq_len(n) %% 3L == 0L),
        Y_2 = as.integer(seq_len(n) %% 5L < 2L),
        A_0 = factor(A_0, levels = c(0, 1)),
        A_1 = factor(A_1, levels = c(0, 1))
    )]
    setcolorder(d, c("id", "W", "L_0", "A_0", "Y_1", "L_1", "A_1", "Y_2"))

    x <- rtmle_init(
        time_grid = 0:2,
        name_id = "id",
        name_outcome = "Y",
        name_competing = NULL,
        name_censoring = NULL,
        minority_threshold = 0
    )
    x$prepared_data <- d
    x$names$name_baseline_covariates <- "W"
    x$names$name_time_covariates <- c("L", "A")
    x$names$name_constant_variables <- NULL
    x
}

make_joint_dynamic_intervention_fixture <- function(repeats = 1L) {
    x <- make_dynamic_intervention_fixture(repeats = repeats)
    x$prepared_data$B_0 <- factor(
        rep(c(0, 1, 0, 1, 0, 0, 1, 1), length.out = NROW(x$prepared_data)),
        levels = c(0, 1)
    )
    x$prepared_data$B_1 <- factor(
        rep(c(0, 0, 1, 0, 0, 1, 0, 1), length.out = NROW(x$prepared_data)),
        levels = c(0, 1)
    )
    data.table::setcolorder(
        x$prepared_data,
        c("id", "W", "L_0", "A_0", "B_0", "Y_1",
          "L_1", "A_1", "B_1", "Y_2")
    )
    x$names$name_time_covariates <- c("L", "A", "B")
    x
}

stop_A_after_L <- function(data, intervention_table, time_node) {
    out <- data.table::copy(data.table::as.data.table(data))
    table_data <- as.data.frame(intervention_table)
    active_rows <- !is.na(table_data$value) & table_data$time_node <= time_node
    for (row in which(active_rows)) {
        action_node <- table_data$time_node[[row]]
        action_variable <- table_data$variable[[row]]
        history_variables <- intersect(
            paste0("L_", 0:max(0, action_node - 1)),
            names(out)
        )
        contraindicated <- if (length(history_variables) == 0L) {
            rep(FALSE, NROW(out))
        } else {
            rowSums(do.call(
                cbind,
                lapply(history_variables, function(v) out[[v]] %in% 1)
            )) > 0
        }
        target_value <- ifelse(
            contraindicated,
            "0",
            as.character(table_data$value[[row]])
        )
        data.table::set(
            out,
            j = action_variable,
            value = factor(target_value, levels = levels(data[[action_variable]]))
        )
    }
    out
}

add_stop_A_regime <- function(
    x,
    name = "A_until_L",
    intervene_function = stop_A_after_L,
    adherence_model_strata = NULL,
    multiple_treatment_factorization = "joint"
) {
    regime(
        x,
        name = name,
        intervention = data.frame(
            time_node = x$intervention_nodes,
            A = factor("1", levels = c("0", "1"))
        ),
        intervene_function = intervene_function,
        multiple_treatment_factorization = multiple_treatment_factorization,
        adherence_model_strata = adherence_model_strata,
        verbose = FALSE
    )
}

add_joint_stop_regime <- function(
    x,
    name = "AB_until_L",
    adherence_model_strata = NULL,
    multiple_treatment_factorization = "joint"
) {
    regime(
        x,
        name = name,
        intervention = data.frame(
            time_node = x$intervention_nodes,
            A = factor("1", levels = c("0", "1")),
            B = factor("0", levels = c("0", "1"))
        ),
        intervene_function = stop_A_after_L,
        multiple_treatment_factorization = multiple_treatment_factorization,
        adherence_model_strata = adherence_model_strata,
        verbose = FALSE
    )
}

formula_rhs_variables <- function(character_formula) {
    all.vars(stats::formula(character_formula)[[3L]])
}

test_that("static interventions preserve numeric treatment coding", {
    data <- data.table(id = 1:2, A_0 = c(0L, 1L))
    intervention_table <- data.table(
        time_node = 0L,
        variable = "A_0",
        value = factor("0", levels = c("0", "1"))
    )
    intervened <- intervene(data, intervention_table, time_node = 0L)
    expect_identical(intervened$A_0, c(0L, 0L))
})

test_that("an intervention function returns data only", {
    old_style_intervention <- function(data, intervention_table, time_node) {
        list(
            data = stop_A_after_L(data, intervention_table, time_node),
            adherence_model_strata = NULL
        )
    }
    expect_error(
        regime(
            make_dynamic_intervention_fixture(),
            name = "Old_style_result",
            intervention = data.frame(
                time_node = 0:1,
                A = factor("1", levels = c("0", "1"))
            ),
            intervene_function = old_style_intervention,
            verbose = FALSE
        ),
        "must return the intervention-updated data"
    )
})

test_that("adherence stratification is dynamic-only", {
    x <- make_dynamic_intervention_fixture()
    expect_error(
        regime(
            x,
            name = "Static_with_strata",
            intervention = data.frame(
                time_node = 0:1,
                A = factor("1", levels = c("0", "1"))
            ),
            adherence_model_strata = "L",
            verbose = FALSE
        ),
        "only for dynamic treatment rules"
    )
})

test_that("the removed dynamic propensity interface is rejected", {
    expect_error(
        regime(
            make_dynamic_intervention_fixture(),
            name = "Old_dynamic_interface",
            intervention = data.frame(
                time_node = 0:1,
                A = factor("1", levels = c("0", "1"))
            ),
            intervene_function = stop_A_after_L,
            dynamic_propensity_instructions = list(strata = "L"),
            verbose = FALSE
        ),
        "removed: dynamic_propensity_instructions.*adherence_model_strata"
    )
})

test_that("dynamic rules transform the response without adding current covariates", {
    x <- make_dynamic_intervention_fixture()
    original_data <- data.table::copy(x$prepared_data)
    x <- add_stop_A_regime(x)
    x <- regime(
        x,
        name = "Static_always_A",
        intervention = data.frame(
            time_node = x$intervention_nodes,
            A = factor("1", levels = c("0", "1"))
        ),
        verbose = FALSE
    )

    expected_match <- matrix(
        c(
            1, 0, 1, 1, 1, 0, 1, 0,
            1, 0, 0, 1, 1, 0, 0, 0
        ),
        ncol = 2
    )
    colnames(expected_match) <- c("A_0", "A_1")
    expect_equal(x$regimes$A_until_L$intervention_match, expected_match)

    x <- model_formula(x, verbose = FALSE)
    dynamic_formula <- x$models$time_1$A_until_L$A_1$formula
    static_formula <- x$models$time_1$Static_always_A$A_1$formula
    expect_match(dynamic_formula, "^\\.rtmle_adherence_A_1 ~")
    expect_true("L_0" %in% formula_rhs_variables(dynamic_formula))
    expect_false("L_1" %in% formula_rhs_variables(dynamic_formula))
    expect_false("L_1" %in% formula_rhs_variables(static_formula))
    expect_identical(x$prepared_data, original_data)
})

test_that("adherence strata resolve to pre-decision history", {
    x <- add_stop_A_regime(
        make_dynamic_intervention_fixture(),
        adherence_model_strata = "L"
    )
    node_one <- rtmle:::evaluate_intervention(
        x$regimes$A_until_L,
        x$prepared_data,
        x$regimes$A_until_L$intervention_table,
        time_node = 1
    )$adherence_model_strata
    expect_identical(node_one, "L_0")
    x <- model_formula(x, Markov = "L", verbose = FALSE)
    node_one_rhs <- formula_rhs_variables(
        x$models$time_1$A_until_L$A_1$formula
    )
    expect_true("L_0" %in% node_one_rhs)
    expect_false("L_1" %in% node_one_rhs)
})

test_that("current-node adherence strata are rejected", {
    x <- make_dynamic_intervention_fixture()
    expect_error(
        add_stop_A_regime(
            x,
            name = "A_current_stratum",
            adherence_model_strata = "L_1"
        ),
        "must use variables known before.*node 0"
    )
})

test_that("strata fit separate adherence models", {
    x <- add_stop_A_regime(
        make_dynamic_intervention_fixture(),
        name = "A_stratified",
        adherence_model_strata = "L"
    )
    x <- model_formula(x, verbose = FALSE)
    fit_sizes <- integer()
    stratified_learner <- function(character_formula, data, intervened_data, ...) {
        fit_sizes <<- c(fit_sizes, NROW(data))
        list(
            predicted_values = rep(0.4, NROW(intervened_data)),
            fit_summary = 0.4
        )
    }
    x <- rtmle:::intervention_probabilities(
        x = x,
        regime_name = "A_stratified",
        max_intervention_node = 1,
        refit = TRUE,
        learner = parse_learners(list(name = "stratified", fun = stratified_learner)),
        seed = 23,
        progressbar = 0,
        save_fitted_objects = FALSE
    )
    expect_length(fit_sizes, 3L)
    expect_equal(
        x$regimes$A_stratified$intervention_probs,
        matrix(c(
            rep(0.4, 8),
            1, 1, rep(0.4, 6)
        ), nrow = 8, ncol = 2,
        dimnames = list(NULL, c("A_0", "A_1")))
    )
})

test_that("multiple treatment factorization is independent of adherence stratification", {
    x <- add_joint_stop_regime(
        make_joint_dynamic_intervention_fixture(),
        name = "AB_sequential",
        multiple_treatment_factorization = "sequential"
    )
    x <- model_formula(x, verbose = FALSE)
    expect_length(x$models$time_0$AB_sequential, 2L)
    expect_true("L_0" %in% formula_rhs_variables(
        x$models$time_1$AB_sequential[[1L]]$A_1$formula
    ))
    expect_true(".rtmle_adherence_A_1" %in% all.vars(
        stats::formula(x$models$time_1$AB_sequential[[1L]]$A_1$formula)
    ))
})

test_that("a dynamic rule without strata uses a pooled adherence model", {
    x <- add_stop_A_regime(
        make_dynamic_intervention_fixture(),
        name = "Invalid_missing_instruction",
        adherence_model_strata = NULL
    )
    x <- model_formula(x, verbose = FALSE)
    x$censoring_use_regime <- "Invalid_missing_instruction"
    observed_adherence <- list()
    learner <- parse_learners(list(name = "pooled", fun = function(
        character_formula, data, intervened_data, ...
    ) {
        observed_adherence[[length(observed_adherence) + 1L]] <<-
            data[["rtmle_outcome"]]
        list(
            predicted_values = rep(0.25, NROW(intervened_data)),
            fit_summary = 0.25
        )
    }))
    x <- rtmle:::intervention_probabilities(
        x = x,
        regime_name = "Invalid_missing_instruction",
        max_intervention_node = 1,
        refit = TRUE,
        learner = learner,
        seed = 13,
        progressbar = 0,
        save_fitted_objects = FALSE
    )
    expect_length(observed_adherence, 2L)
    expect_equal(observed_adherence[[2L]], c(1L, 1L, 0L, 1L, 1L, 1L, 0L, 0L))
})

test_that("dynamic adherence probabilities are fitted in G", {
    x <- add_stop_A_regime(make_dynamic_intervention_fixture())
    x <- model_formula(x, verbose = FALSE)
    observed_adherence <- list()
    adherence_learner <- function(character_formula, data, intervened_data, ...) {
        observed_adherence[[length(observed_adherence) + 1L]] <<-
            data[["rtmle_outcome"]]
        list(
            predicted_values = rep(0.25, NROW(intervened_data)),
            fit_summary = 0.25
        )
    }
    x <- rtmle:::intervention_probabilities(
        x = x,
        regime_name = "A_until_L",
        max_intervention_node = 1,
        refit = TRUE,
        learner = parse_learners(list(name = "adherence", fun = adherence_learner)),
        seed = 19,
        progressbar = 0,
        save_fitted_objects = FALSE
    )
    expect_equal(
        observed_adherence[[2L]],
        c(1L, 1L, 0L, 1L, 1L, 1L, 0L, 0L)
    )
    expect_equal(
        x$regimes$A_until_L$intervention_probs,
        matrix(0.25, nrow = 8, ncol = 2,
               dimnames = list(NULL, c("A_0", "A_1")))
    )
})

test_that("run_rtmle consumes dynamic adherence models", {
    x <- add_stop_A_regime(
        make_dynamic_intervention_fixture(repeats = 20L),
        name = "A_until_L"
    )
    x <- target(
        x,
        name = "Outcome_risk",
        estimator = "tmle",
        regimes = "A_until_L"
    )
    x <- model_formula(x, verbose = FALSE)
    x <- suppressWarnings(run_rtmle(
        x,
        learner = "learn_glm",
        time_horizon = 2,
        refit = TRUE,
        verbose = FALSE
    ))
    expect_equal(NROW(x$estimate$Main_analysis), 1L)
    expect_true(is.finite(x$estimate$Main_analysis$Estimate))
    expect_gte(x$estimate$Main_analysis$Estimate, 0)
    expect_lte(x$estimate$Main_analysis$Estimate, 1)
})
