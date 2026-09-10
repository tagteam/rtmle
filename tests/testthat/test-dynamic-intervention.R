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


stop_A_after_L <- function(data, intervention_table, time_node) {
    intervened_data <- data.table::copy(data.table::as.data.table(data))
    table_data <- as.data.frame(intervention_table)
    active_rows <- !is.na(table_data$value) & table_data$time_node <= time_node

    for (row in which(active_rows)) {
        action_node <- table_data$time_node[[row]]
        action_variable <- table_data$variable[[row]]
        history_variables <- intersect(
            paste0("L_", 0:action_node),
            names(intervened_data)
        )
        contraindicated <- Reduce(
            `|`,
            lapply(history_variables, function(v) intervened_data[[v]] %in% 1)
        )
        target_value <- ifelse(
            contraindicated,
            "0",
            as.character(table_data$value[[row]])
        )
        data.table::set(
            intervened_data,
            j = action_variable,
            value = factor(target_value, levels = levels(data[[action_variable]]))
        )
    }

    intervened_data
}


fixed_stop_A_instructions <- function(data, intervention_table, time_node) {
    history_variables <- intersect(
        paste0("L_", 0:time_node),
        names(data)
    )
    contraindicated <- Reduce(
        `|`,
        lapply(history_variables, function(v) data[[v]] %in% 1)
    )
    action <- paste0("A_", time_node)
    stats::setNames(
        list(list(
            mode = "fixed",
            probability = ifelse(contraindicated, 1, NA_real_)
        )),
        action
    )
}


adherence_stop_A_instructions <- function(data, intervention_table, time_node,
                                          stratify_by = NULL) {
    action <- paste0("A_", time_node)
    stats::setNames(
        list(list(mode = "adherence", stratify_by = stratify_by)),
        action
    )
}


stop_A_propensity_variables <- function(data, intervention_table, time_node) {
    action <- paste0("A_", time_node)
    # L is cumulative in this fixture, so its current value summarizes its
    # history. Naming L_time_node also opts into the ordering L_k -> A_k.
    stats::setNames(list(paste0("L_", time_node)), action)
}


add_stop_A_protocol <- function(x,
                                name = "A_until_L",
                                intervene_function = stop_A_after_L,
                                propensity_instructions = fixed_stop_A_instructions,
                                propensity_variables = stop_A_propensity_variables) {
    protocol(
        x,
        name = name,
        intervention = data.frame(
            time_node = x$intervention_nodes,
            A = factor("1", levels = c("0", "1"))
        ),
        intervene_function = intervene_function,
        propensity_instructions = propensity_instructions,
        propensity_variables = propensity_variables,
        verbose = FALSE
    )
}


formula_rhs_variables <- function(character_formula) {
    all.vars(stats::formula(character_formula)[[3L]])
}


test_that("protocol propensity metadata are parsed and validated", {
    expect_false("dynamic_intervention_instructions" %in% getNamespaceExports("rtmle"))
    expect_false(exists("intervention_result", envir = asNamespace("rtmle"),
                        inherits = FALSE))
    old_style_intervention <- function(data, intervention_table, time_node) {
        list(
            data = stop_A_after_L(data, intervention_table, time_node),
            propensity_instructions = NULL,
            propensity_variables = NULL
        )
    }
    expect_error(
        protocol(
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
    x <- make_dynamic_intervention_fixture()
    x <- protocol(
        x,
        name = "A_until_L_metadata",
        intervention = data.frame(
            time_node = x$intervention_nodes,
            A = factor("1", levels = c("0", "1"))
        ),
        intervene_function = stop_A_after_L,
        propensity_instructions = list(
            A_0 = list(mode = "fixed", probability = rep(NA_real_, 8))
        ),
        propensity_variables = list(A_0 = "L_0"),
        verbose = FALSE
    )

    expect_equal(
        x$protocols$A_until_L_metadata$propensity_instructions$A_0$mode,
        "fixed"
    )
    expect_equal(
        x$protocols$A_until_L_metadata$propensity_instructions$A_0$probability,
        rep(NA_real_, 8)
    )
    expect_identical(
        x$protocols$A_until_L_metadata$propensity_variables,
        list(A_0 = "L_0")
    )

    expect_error(protocol(
        make_dynamic_intervention_fixture(),
        name = "Invalid_probability_length",
        intervention = data.frame(
            time_node = 0:1,
            A = factor("1", levels = c("0", "1"))
        ),
        intervene_function = stop_A_after_L,
        propensity_instructions = list(
            A_0 = list(mode = "fixed", probability = c(NA_real_, 1))
        ),
        verbose = FALSE
    ), "length")
    expect_error(protocol(
        make_dynamic_intervention_fixture(),
        name = "Invalid_probability_range",
        intervention = data.frame(
            time_node = 0:1,
            A = factor("1", levels = c("0", "1"))
        ),
        intervene_function = stop_A_after_L,
        propensity_instructions = list(
            A_0 = list(mode = "fixed", probability = c(NA_real_, 1.01, rep(NA_real_, 6)))
        ),
        verbose = FALSE
    ), "\\[0, 1\\]")
    expect_error(protocol(
        make_dynamic_intervention_fixture(),
        name = "Invalid_probability_nan",
        intervention = data.frame(
            time_node = 0:1,
            A = factor("1", levels = c("0", "1"))
        ),
        intervene_function = stop_A_after_L,
        propensity_instructions = list(
            A_0 = list(mode = "fixed", probability = c(NA_real_, NaN, rep(NA_real_, 6)))
        ),
        verbose = FALSE
    ), "NaN")
})


test_that("intervention metadata are validated against the supplied data", {
    unknown_propensity_variable <- function(data, intervention_table, time_node) {
        "not_in_data"
    }
    x <- make_dynamic_intervention_fixture()

    expect_error(
        add_stop_A_protocol(
            x,
            name = "Invalid_rule",
            propensity_variables = unknown_propensity_variable
        ),
        "Unknown.*propensity_variables"
    )
})


test_that("a contraindication rule controls formulas and row-wise adherence", {
    x <- make_dynamic_intervention_fixture()
    original_data <- data.table::copy(x$prepared_data)
    x <- add_stop_A_protocol(x)
    x <- protocol(
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
            1, 0, 1, 0, 1, 0, 0, 0
        ),
        ncol = 2
    )
    colnames(expected_match) <- c("A_0", "A_1")
    expect_equal(x$protocols$A_until_L$intervention_match, expected_match)

    # In particular, stopping A after L is adherent while the rare continuation
    # of A after L is not. This must be a row-wise comparison, not `%in%`.
    expect_equal(x$protocols$A_until_L$intervention_match[c(3, 4), "A_1"], c(1, 0))

    x <- model_formula(x, verbose = FALSE)
    dynamic_rhs <- formula_rhs_variables(
        x$models$time_1$A_until_L$A_1$formula
    )
    static_rhs <- formula_rhs_variables(
        x$models$time_1$Static_always_A$A_1$formula
    )
    expect_true("L_1" %in% dynamic_rhs)
    expect_false("L_1" %in% static_rhs)
    expect_identical(x$prepared_data, original_data)
})


test_that("a current-node propensity variable replaces its Markov history", {
    x <- make_dynamic_intervention_fixture()
    x <- add_stop_A_protocol(x)
    x <- model_formula(x, Markov = "L", verbose = FALSE)

    rhs <- formula_rhs_variables(x$models$time_1$A_until_L$A_1$formula)
    expect_true("L_1" %in% rhs)
    expect_false("L_0" %in% rhs)
})


test_that("fixed propensity instructions skip their deterministic rows", {
    x <- make_dynamic_intervention_fixture()
    x <- add_stop_A_protocol(x)
    x <- model_formula(x, verbose = FALSE)
    x$censoring_use_protocol <- "A_until_L"

    fit_sizes <- integer()
    quarter_learner <- function(character_formula, data, intervened_data, ...) {
        fit_sizes <<- c(fit_sizes, NROW(data))
        list(
            predicted_values = rep(0.25, NROW(intervened_data)),
            fit_summary = 0.25
        )
    }
    learner <- parse_learners(list(name = "quarter", fun = quarter_learner))

    x <- rtmle:::intervention_probabilities(
        x = x,
        protocol_name = "A_until_L",
        max_intervention_node = 1,
        refit = TRUE,
        learner = learner,
        seed = 11,
        progressbar = 0,
        save_fitted_objects = FALSE
    )

    expected_probabilities <- cbind(
        A_0 = c(1, 1, rep(0.25, 6)),
        A_1 = c(rep(1, 4), rep(0.25, 4))
    )
    expect_equal(x$protocols$A_until_L$intervention_probs, expected_probabilities)
    expect_equal(fit_sizes, c(6L, 4L))
    expect_equal(
        x$protocols$A_until_L$cumulative_intervention_probs[, 2],
        c(1, 1, 0.25, 0.25, rep(0.0625, 4))
    )
})


test_that("static fixed probabilities are subset for at-risk fits", {
    x <- make_dynamic_intervention_fixture()
    x$followup <- data.table::data.table(
        last_interval = c(1, 0, 1, 0, 1, 1, 0, 1)
    )
    x <- add_stop_A_protocol(
        x,
        name = "A_until_L_static",
        propensity_instructions = list(
            A_0 = list(
                mode = "fixed",
                probability = c(1, 1, rep(NA_real_, 6))
            ),
            A_1 = list(
                mode = "fixed",
                probability = c(1, 1, 1, 1, rep(NA_real_, 4))
            )
        ),
        propensity_variables = NULL
    )
    x <- model_formula(x, verbose = FALSE)
    learner <- parse_learners(list(
        name = "quarter",
        fun = function(character_formula, data, intervened_data, ...) {
            list(
                predicted_values = rep(0.25, NROW(intervened_data)),
                fit_summary = 0.25
            )
        }
    ))
    x <- rtmle:::intervention_probabilities(
        x = x,
        protocol_name = "A_until_L_static",
        max_intervention_node = 1,
        refit = TRUE,
        learner = learner,
        seed = 17,
        progressbar = 0,
        save_fitted_objects = FALSE
    )

    expect_equal(
        x$protocols$A_until_L_static$intervention_probs,
        rbind(
            c(A_0 = 1, A_1 = 1),
            c(A_0 = 1, A_1 = NA),
            c(A_0 = 0.25, A_1 = 1),
            c(A_0 = 0.25, A_1 = NA),
            c(A_0 = 0.25, A_1 = 0.25),
            c(A_0 = 0.25, A_1 = 0.25),
            c(A_0 = 0.25, A_1 = NA),
            c(A_0 = 0.25, A_1 = 0.25)
        )
    )
})


test_that("changing a nominal treatment value requires a propensity instruction", {
    stop_without_instruction <- function(data, intervention_table, time_node) {
        stop_A_after_L(data, intervention_table, time_node)
    }
    x <- make_dynamic_intervention_fixture()
    x <- add_stop_A_protocol(
        x,
        name = "Invalid_missing_instruction",
        intervene_function = stop_without_instruction,
        propensity_instructions = NULL,
        propensity_variables = NULL
    )
    x <- model_formula(x, verbose = FALSE)
    x$censoring_use_protocol <- "Invalid_missing_instruction"
    unused_learner <- function(...) {
        stop("validation must happen before fitting")
    }
    learner <- parse_learners(list(name = "unused", fun = unused_learner))

    expect_error(
        rtmle:::intervention_probabilities(
            x = x,
            protocol_name = "Invalid_missing_instruction",
            max_intervention_node = 1,
            refit = TRUE,
            learner = learner,
            seed = 13,
            progressbar = 0,
            save_fitted_objects = FALSE
        ),
        "changes a nominal treatment value.*without a propensity instruction"
    )
})


test_that("adherence instructions fit the dynamic match outcome", {
    x <- make_dynamic_intervention_fixture()
    x <- add_stop_A_protocol(
        x,
        name = "A_until_L_adherence",
        propensity_instructions = adherence_stop_A_instructions
    )
    x <- model_formula(x, verbose = FALSE)
    adherence_formula <- x$models$time_1$A_until_L_adherence$A_1$formula
    expect_match(adherence_formula, "^\\.rtmle_adherence_A_1 ~")
    expect_match(adherence_formula, "L_1", fixed = TRUE)

    observed_adherence <- list()
    fit_sizes <- integer()
    adherence_learner <- function(character_formula, data, intervened_data, ...) {
        fit_sizes <<- c(fit_sizes, NROW(data))
        observed_adherence[[length(observed_adherence) + 1L]] <<-
            data[["rtmle_outcome"]]
        list(
            predicted_values = rep(0.25, NROW(intervened_data)),
            fit_summary = 0.25
        )
    }
    x <- rtmle:::intervention_probabilities(
        x = x,
        protocol_name = "A_until_L_adherence",
        max_intervention_node = 1,
        refit = TRUE,
        learner = parse_learners(list(name = "adherence", fun = adherence_learner)),
        seed = 19,
        progressbar = 0,
        save_fitted_objects = FALSE
    )

    expect_equal(fit_sizes, c(8L, 8L))
    expect_equal(
        observed_adherence[[2L]],
        c(1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L)
    )
    expect_equal(
        x$protocols$A_until_L_adherence$intervention_probs,
        matrix(0.25, nrow = 8, ncol = 2,
               dimnames = list(NULL, c("A_0", "A_1")))
    )
})


test_that("adherence instructions can fit separate contraindication strata", {
    stratified_rule <- function(data, intervention_table, time_node) {
        stop_A_after_L(data, intervention_table, time_node)
    }
    stratified_instructions <- function(data, intervention_table, time_node) {
        adherence_stop_A_instructions(
            data,
            intervention_table,
            time_node,
            stratify_by = paste0("L_", time_node)
        )
    }
    x <- make_dynamic_intervention_fixture()
    x <- add_stop_A_protocol(
        x,
        name = "A_until_L_stratified",
        intervene_function = stratified_rule,
        propensity_instructions = stratified_instructions
    )
    x <- model_formula(x, verbose = FALSE)
    expect_match(
        x$models$time_1$A_until_L_stratified$A_1$formula,
        "L_1",
        fixed = TRUE
    )

    stratified_learner <- function(character_formula, data, intervened_data, ...) {
        list(
            predicted_values = rep(mean(data[["rtmle_outcome"]]),
                                   NROW(intervened_data)),
            fit_summary = mean(data[["rtmle_outcome"]])
        )
    }
    x <- rtmle:::intervention_probabilities(
        x = x,
        protocol_name = "A_until_L_stratified",
        max_intervention_node = 1,
        refit = TRUE,
        learner = parse_learners(list(name = "stratified", fun = stratified_learner)),
        seed = 23,
        progressbar = 0,
        save_fitted_objects = FALSE
    )

    expect_equal(
        x$protocols$A_until_L_stratified$intervention_probs[, "A_1"],
        c(rep(0.75, 4), rep(0.5, 4))
    )
})


test_that("fully deterministic treatment nodes do not call the learner", {
    x <- make_dynamic_intervention_fixture()
    x$prepared_data[, `:=`(L_0 = 1, L_1 = 1)]
    x <- add_stop_A_protocol(x)
    x <- model_formula(x, verbose = FALSE)
    x$censoring_use_protocol <- "A_until_L"

    unused_learner <- function(...) {
        stop("the learner must not be called for a fully overridden node")
    }
    learner <- parse_learners(list(name = "unused", fun = unused_learner))

    x <- rtmle:::intervention_probabilities(
        x = x,
        protocol_name = "A_until_L",
        max_intervention_node = 1,
        refit = TRUE,
        learner = learner,
        seed = 12,
        progressbar = 0,
        save_fitted_objects = FALSE
    )

    expect_true(all(x$protocols$A_until_L$intervention_probs == 1))
    expect_true(all(x$protocols$A_until_L$cumulative_intervention_probs == 1))
})


test_that("re-registering a protocol invalidates its derived probability caches", {
    x <- make_dynamic_intervention_fixture()
    x <- add_stop_A_protocol(x)
    old_match <- x$protocols$A_until_L$intervention_match
    x$protocols$A_until_L$intervention_probs <- matrix(0.5, nrow = 8, ncol = 2)
    x$protocols$A_until_L$cumulative_intervention_probs <- matrix(0.25, nrow = 8, ncol = 2)
    x$protocols$A_until_L$ipw_last_nodes <- c(node_0 = "A_0", node_1 = "A_1")
    x$protocols$A_until_L$intervention_last_nodes <- c(node_0 = "A_0", node_1 = "A_1")

    x <- protocol(
        x,
        name = "A_until_L",
        intervention = data.frame(
            time_node = x$intervention_nodes,
            A = factor("0", levels = c("0", "1"))
        ),
        verbose = FALSE
    )

    expect_null(x$protocols$A_until_L$intervention_probs)
    expect_null(x$protocols$A_until_L$cumulative_intervention_probs)
    expect_null(x$protocols$A_until_L$ipw_last_nodes)
    expect_null(x$protocols$A_until_L$intervention_last_nodes)
    expect_false(identical(x$protocols$A_until_L$intervention_match, old_match))

    expected_match <- matrix(
        c(
            1, 0, 0, 0, 0, 1, 0, 1,
            1, 0, 0, 0, 0, 0, 0, 1
        ),
        ncol = 2
    )
    colnames(expected_match) <- c("A_0", "A_1")
    expect_equal(x$protocols$A_until_L$intervention_match, expected_match)
})


test_that("data-only intervention functions remain supported", {
    x <- make_dynamic_intervention_fixture()
    x <- add_stop_A_protocol(
        x,
        name = "Data_only_A_until_L",
        intervene_function = stop_A_after_L,
        propensity_instructions = NULL,
        propensity_variables = NULL
    )

    expected_match <- matrix(
        c(
            1, 0, 1, 1, 1, 0, 1, 0,
            1, 0, 1, 0, 1, 0, 0, 0
        ),
        ncol = 2
    )
    colnames(expected_match) <- c("A_0", "A_1")
    expect_equal(x$protocols$Data_only_A_until_L$intervention_match, expected_match)

    x <- model_formula(x, verbose = FALSE)
    expect_false(
        "L_1" %in% formula_rhs_variables(
            x$models$time_1$Data_only_A_until_L$A_1$formula
        )
    )
})


test_that("run_rtmle consumes protocol propensity instructions in G and Q", {
    x <- make_dynamic_intervention_fixture(repeats = 20L)
    x <- add_stop_A_protocol(x)
    x <- target(
        x,
        name = "Outcome_risk",
        estimator = "tmle",
        protocols = "A_until_L"
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
    expect_true(all(
        x$protocols$A_until_L$intervention_probs[x$prepared_data$L_1 %in% 1, "A_1"] == 1
    ))
})


test_that("run_rtmle fits pooled adherence probabilities", {
    x <- make_dynamic_intervention_fixture(repeats = 20L)
    x <- add_stop_A_protocol(
        x,
        name = "A_until_L_adherence",
        propensity_instructions = adherence_stop_A_instructions
    )
    x <- target(
        x,
        name = "Outcome_risk",
        estimator = "tmle",
        protocols = "A_until_L_adherence"
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
    expect_true(any(
        x$protocols$A_until_L_adherence$intervention_probs[, "A_1"] < 1
    ))
})


test_that("simulated_cohort supports stopping A after bleeding", {
    data(simulated_cohort, package = "rtmle")
    ld <- register_format(simulated_cohort)
    x <- rtmle_init(
        time_grid = seq(0, 20, 4),
        name_id = "id",
        name_outcome = "stroke",
        name_competing = "death",
        name_censoring = NULL
    )
    x <- add_long_data(
        x,
        outcome_data = ld$timevar_data$stroke[!duplicated(id)],
        censored_data = NULL,
        competing_data = ld$timevar_data$death,
        timevar_data = ld$timevar_data[c("bleeding", "changeSBP", "A", "B")]
    )
    x <- add_baseline_data(x, data = ld$baseline_data)
    x <- long_to_wide(x, start_followup_date = 0)
    x <- prepare_rtmle_data(x)
    original_data <- data.table::copy(x$prepared_data)

    has_bled_by <- function(data, node) {
            history <- intersect(paste0("bleeding_", 0:node), names(data))
            if (length(history) == 0) {
                return(rep(FALSE, NROW(data)))
            }
            bleeding_history <- do.call(
                cbind,
                lapply(history, function(v) data[[v]] %in% c(1, "1"))
            )
            rowSums(bleeding_history) > 0
    }

    stop_A_after_bleeding <- function(data, intervention_table, time_node) {
        intervened_data <- intervene(data, intervention_table, time_node)
        for (node in unique(intervention_table$time_node)) {
            action <- paste0("A_", node)
            intervened_data[[action]][has_bled_by(data, node)] <- "0"
        }
        intervened_data
    }

    bleeding_propensity_instructions <- function(data, intervention_table,
                                                  time_node) {
        current_action <- paste0("A_", time_node)
        stats::setNames(
            list(list(
                mode = "fixed",
                probability = ifelse(
                    has_bled_by(data, time_node),
                    1,
                    NA_real_
                )
            )),
            current_action
        )
    }

    bleeding_propensity_variables <- function(data, intervention_table,
                                              time_node) {
        current_action <- paste0("A_", time_node)
        stats::setNames(
            list(paste0("bleeding_", time_node)),
            current_action
        )
    }

    x <- protocol(
        x,
        name = "A_until_bleeding",
        intervention = data.frame(
            time_node = x$intervention_nodes,
            A = factor("1", levels = c("0", "1"))
        ),
        intervene_function = stop_A_after_bleeding,
        propensity_instructions = bleeding_propensity_instructions,
        propensity_variables = bleeding_propensity_variables,
        verbose = FALSE
    )

    # Independently reconstruct the expected dynamic adherence matrix. A
    # continuation of A after bleeding is not adherent, whereas stopping is.
    expected_match <- matrix(
        0L,
        nrow = NROW(x$prepared_data),
        ncol = 2L,
        dimnames = list(NULL, c("A_0", "A_1"))
    )
    previous <- rep(TRUE, NROW(x$prepared_data))
    bled <- rep(FALSE, NROW(x$prepared_data))
    for (node in 0:1) {
        current_bleeding <- x$prepared_data[[paste0("bleeding_", node)]] %in% 1
        current_bleeding[is.na(current_bleeding)] <- FALSE
        bled <- bled | current_bleeding
        desired <- ifelse(bled, 0, 1)
        observed <- x$prepared_data[[paste0("A_", node)]]
        same <- !is.na(observed) & observed == desired
        previous <- previous & same
        expected_match[, node + 1L] <- as.integer(previous)
    }
    expect_equal(
        x$protocols$A_until_bleeding$intervention_match[, c("A_0", "A_1")],
        expected_match
    )
    at_risk_1 <- x$followup$last_interval >= 1
    bleeding_1 <- !is.na(x$prepared_data$bleeding_1) &
        x$prepared_data$bleeding_1 %in% 1
    expect_gt(sum(at_risk_1 & bleeding_1), 0)
    expect_true(all(
        x$protocols$A_until_bleeding$intervention_match[
            at_risk_1 & bleeding_1 & x$prepared_data$A_0 %in% 1 &
                x$prepared_data$A_1 %in% 0,
            "A_1"
        ] == 1
    ))
    expect_false(any(
        x$protocols$A_until_bleeding$intervention_match[
            at_risk_1 & bleeding_1 & x$prepared_data$A_1 %in% 1,
            "A_1"
        ] == 1
    ))

    x <- model_formula(x, verbose = FALSE)
    dynamic_formula <- x$models$time_1$A_until_bleeding$A_1$formula
    expect_true(grepl("bleeding_1", dynamic_formula, fixed = TRUE))
    expect_false(grepl("bleeding_2", dynamic_formula, fixed = TRUE))
    expect_identical(x$prepared_data, original_data)

    # A sentinel learner makes the deterministic instruction visible: bleeding
    # rows receive one and only non-bleeding rows are passed to the learner.
    fit_sizes <- integer()
    sentinel_learner <- function(character_formula, data, intervened_data, ...) {
        fit_sizes <<- c(fit_sizes, NROW(data))
        list(
            predicted_values = rep(0.4, NROW(intervened_data)),
            fit_summary = 0.4
        )
    }
    x <- rtmle:::intervention_probabilities(
        x = x,
        protocol_name = "A_until_bleeding",
        max_intervention_node = 1,
        refit = TRUE,
        learner = parse_learners(list(name = "sentinel", fun = sentinel_learner)),
        seed = 104,
        progressbar = 0,
        save_fitted_objects = FALSE
    )
    intervention_probs <- x$protocols$A_until_bleeding$intervention_probs
    expect_true(all(intervention_probs[at_risk_1 & bleeding_1, "A_1"] == 1))
    expect_true(any(intervention_probs[at_risk_1 & !bleeding_1, "A_1"] == 0.4))
    expect_lt(fit_sizes[[2]], fit_sizes[[1]])

    x <- target(
        x,
        name = "Stroke_risk",
        estimator = "tmle",
        protocols = "A_until_bleeding"
    )
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
