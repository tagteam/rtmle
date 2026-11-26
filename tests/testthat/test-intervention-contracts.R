test_that("regimes reject invalid intervention schedules and arguments", {
    x <- rtmle_init(time_grid = 0:2, name_id = "id", name_outcome = "Y")
    register <- function(...) regime(x, "rule", 1, treatment_variables = "A", ...)
    schedule <- function(table) regime(x, "rule", table)
    action <- factor(1, levels = 0:1)

    expect_error(register(dynamic_propensity_instructions = "joint"), "former interface")
    expect_error(register(unknown = TRUE), "Unused argument.*unknown")
    expect_error(register(multiple_treatment_factorization = "marginal"), "factorization")
    expect_error(register(adherence_model_strata = NA_character_), "character vector")
    expect_error(register(adherence_model_strata = "age"), "only for dynamic")
    expect_error(register(intervene_function = FALSE), "function or the name")
    expect_error(register(intervene_function_args = 3), "named list")
    expect_error(register(intervene_function_args = list(time_node = 1)), "non-reserved")
    expect_error(register(intervene_function_args = list(cutoff = 4)), "requires")
    expect_error(regime(x, "rule", 2, treatment_variables = "A"), "vector of 0s and 1s")
    expect_error(schedule(data.frame(A = action)), "intervention_node")
    expect_error(schedule(data.frame(time_node = 0, A_0 = action)), "without time suffix")
    expect_error(schedule(data.frame(time_node = 2, A = action)), "nodes in the object: 2")
    expect_error(schedule(data.frame(time_node = 0, A = 1)), "must be factors")
    expect_error(schedule(data.frame(time_node = 0, A = factor(1))), "exactly 2 levels")

    # When the node index is explicit, a treatment column may itself be called node.
    dual <- schedule(data.frame(time_node = 0, node = action))
    expect_identical(dual$regimes$rule$treatment_variables, "node")

    x <- register()
    x <- regime(x, "rule_B", treatment_variables = "B", intervention = 0)
    expect_named(x$names$treatment_options, c("A", "B"))
    x$prepared_data <- data.table::data.table(id = 1L)
    expect_error(intervention_match(x, "rule"), "not in x\\$prepared_data: A_0")
})

test_that("direct interventions apply assignments only through the requested node", {
    history <- data.table::data.table(A_0 = 0:1, A_1 = 0:1, bleeding_1 = c(FALSE, TRUE))
    assignments <- data.table::data.table(time_node = 0:1, variable = c("A_0", "A_1"),
                                          value = factor(c("1", "1"), levels = 0:1))
    baseline <- intervene(history, assignments, time_node = 0)
    expect_identical(baseline$A_0, c(1L, 1L))
    expect_identical(baseline$A_1, 0:1)
    next_node <- intervene(history, assignments, time_node = 1)
    expect_identical(next_node$A_0, c(1L, 1L))
    expect_identical(next_node$A_1, c(1L, 1L))

    baseline_rule <- treat_unless(history, assignments, time_node = 0,
                                  contra_indication = "bleeding", action = 0)
    expect_identical(baseline_rule$A_1, 0:1)
})

test_that("intervention callbacks receive the current node and retain the data layout", {
    history <- data.table::data.table(id = 1:2, A_0 = 0:1, L_0 = 0:1)
    assignment <- data.table::data.table(time_node = 0L, variable = "A_0", value = 1L)
    evaluate <- function(callback, strata = NULL, node = 1L, ...) {
        rtmle:::evaluate_intervention(
            list(intervene_function = callback, adherence_model_strata = strata),
            history, assignment, time_node = node, ...
        )
    }
    identity_rule <- function(data, intervention_table) data
    via_dots <- function(data, ...) {
        data$A_0 <- rep(list(...)$time_node, nrow(data))
        data
    }
    via_third_argument <- function(data, intervention_table, decision) {
        data$A_0 <- rep(decision, nrow(data))
        data
    }
    expect_equal(evaluate(via_dots)$data$A_0, c(1L, 1L))
    expect_equal(evaluate(via_third_argument)$data$A_0, c(1L, 1L))
    expect_identical(history$A_0, 0:1)
    expect_equal(evaluate(identity_rule, c("L_0", "id"))$adherence_model_strata,
                 c("L_0", "id"))
    expect_error(evaluate(identity_rule, "L_1"), "not pre-decision")
    expect_error(evaluate(identity_rule, "age"), "Unknown variable.*age")
    expect_error(evaluate(identity_rule, "L", node = Inf), "one finite value")
    expect_error(evaluate(identity_rule, full_n = -1), "one non-negative integer")
    expect_error(evaluate(identity_rule, row_indices = c(1, 3)), "full cohort")
    expect_error(evaluate(NULL), "no intervene_function")
    expect_error(evaluate("missing_treatment_callback"), "missing_treatment_callback")
    expect_error(evaluate(2), "function or the name")
    expect_error(evaluate(function(data, ...) TRUE), "data frame or matrix")
    expect_error(evaluate(function(data, ...) data[1L]), "number and order of rows")
    expect_error(evaluate(function(data, ...) data[, -"A_0"]), "removed required.*A_0")
    expect_error(rtmle:::evaluate_intervention(
        list(intervene_function = identity_rule, intervene_function_args = list(1L)),
        history, assignment, time_node = 0L), "named list")
})

test_that("contraindication rules validate their configuration and preserve treatment types", {
    assignment <- data.table::data.table(variable = "A_0", value = factor(1, levels = 0:1))
    history <- data.table::data.table(A_0 = 0:1, bleeding = c(TRUE, FALSE))
    rule <- function(action = 0, event = "bleeding") {
        treat_unless(history, assignment, 0, event, action)
    }
    expect_error(treat_unless(contra_indication = "bleeding"), "needs both")
    expect_error(treat_unless(history, assignment, contra_indication = "bleeding", action = 0),
                 "must be supplied")
    expect_error(treat_unless(history, assignment, 0, action = 0), "contra_indication and action must be supplied")
    expect_error(rule(event = "bleeding_0"), "without time suffixes")
    expect_error(rule(event = "stroke"), "not present in data: stroke")
    expect_error(rule(action = list(0)), "Argument action must contain")
    expect_error(rule(action = c(B = 0)), "must name each treatment variable")
    expect_error(rule(action = c(0, 1)), "one value or one value for each")
    expect_identical(rule(event = NULL)$A_0, c(1L, 1L))
    expect_identical(rule(action = c(A = 0L))$A_0, c(0L, 1L))

    future_bleeding <- data.table::data.table(A_0 = 0:1, bleeding_1 = c(TRUE, FALSE))
    at_baseline <- treat_unless(future_bleeding, assignment, 0, "bleeding", 0)
    expect_identical(at_baseline$A_0, c(1L, 1L))

    history$A_0 <- c(TRUE, FALSE)
    expect_identical(rule()$A_0, c(FALSE, TRUE))
    history$A_0 <- c(0, 0)
    expect_identical(rule()$A_0, c(0, 1))
    history$bleeding <- c("true", "false")
    expect_identical(rule()$A_0, c(0, 1))
    assignment$value <- 1L
    history$A_0 <- c(0L, 0L)
    expect_identical(rule(action = 0.0)$A_0, c(0L, 1L))
    history$A_0 <- c("0", "1")
    assignment$value <- "1"
    expect_identical(rule(action = "0")$A_0, c("0", "1"))
    assignment$value <- TRUE
    history$A_0 <- c(FALSE, FALSE)
    expect_identical(rule(action = "false")$A_0, c(FALSE, TRUE))
})

test_that("adherence reports missing inputs and handles regimes with no initiators", {
    x <- list(followup = data.frame(id = 1:2, last_interval = c(1L, 1L)),
              names = list(outcome = "Y", competing = "D", censoring = "C"))
    expect_error(adherence(x), "No regime")
    x$regimes <- list(empty = list())
    expect_error(adherence(x), "No regime")
    x$regimes$untreated <- list(intervention_match = matrix(0, nrow = 2, ncol = 2))
    expect_error(adherence(x, ""), "non-empty character vector")
    expect_error(adherence(x, "unknown"), "does not contain the prepared data")
    expect_error(summary_adherence(x), "requires x\\$prepared_data")
    empty <- adherence(x)
    expect_named(empty, c("id", "regime", "time_nonadherence", "event_nonadherence"))
    expect_equal(nrow(empty), 0)
    expect_error(summary(empty), "has no summary")

    x$prepared_data <- data.frame(id = 1:2, Y_1 = 0:1)
    x$names$competing <- NULL
    expect_error(adherence(x, "unknown"), "no element intervention_match")
    counts <- summary_adherence(x)
    expect_identical(counts$n_initiated, c(0L, 0L))
    expect_equal(counts$n_outcome, c(0L, 0L))
    expect_equal(counts$n_death, c(0L, 0L))
    x$regimes$untreated$intervention_match <- matrix(integer(), nrow = 2, ncol = 0)
    expect_equal(nrow(adherence(x)), 0)
    expect_equal(nrow(summary_adherence(x)), 0)

    x$names$id <- "pid"
    expect_error(adherence(x), "identifier column pid")
    x$followup$last_interval <- NULL
    expect_error(adherence(x), "last_interval")
    expect_error(summary_adherence(x), "last_interval")
})
