continuous_intervention_policy <- function(data, regime, node) {
    if (!nrow(data)) return(list(value = numeric(), natural = logical()))
    if (!isTRUE(regime$dynamic_intervention)) return(list(
        value = rep(regime$static_intervention, nrow(data)), natural = rep(FALSE, nrow(data))))

    columns <- names(data)
    indexed <- grepl("_[0-9]+$", columns)
    index <- rep(NA_integer_, length(columns))
    index[indexed] <- as.integer(sub("^.*_", "", columns[indexed]))
    internal <- grepl("^(propensity|survival_censoring|time_pooled|event_pooled)_", columns)
    # Fixed baseline predictors may themselves end in digits (e.g. year_2024).
    # They remain available at every node rather than being mistaken for future
    # event-history columns.
    baseline <- regime$continuous_baseline_covariates
    history <- data.table::copy(data[, columns[!internal &
        (!indexed | index <= node | columns %in% baseline)], with = FALSE])
    # Supply the discrete intervention evaluator with the treatment nodes it
    # expects, while restricting its data to the history available at this node.
    nodes <- 0:node
    table <- data.table::data.table(time_node = nodes, variable = paste0("A_", nodes),
                                    value = factor(rep(regime$static_intervention, length(nodes)), levels = c(0, 1)))
    column <- paste0("A_", node)
    observed <- as.numeric(as.character(history[[column]]))
    evaluate <- function(current) {
        input <- data.table::copy(history)
        data.table::set(input, j = column, value = rep(current, nrow(input)))
        # The discrete evaluator's strata resolver uses discrete node timing.
        # Continuous strata are resolved separately when fitting adherence.
        rule <- regime
        rule$adherence_model_strata <- NULL
        intervened <- evaluate_intervention(rule, input, table, node)$data
        if (!identical(intervened$id, input$id))
            stop("The intervention function must preserve subject identifiers and row order.")
        values <- as.character(intervened[[column]])
        if (length(values) != nrow(input) || anyNA(values) || !all(values %in% c("0", "1")))
            stop("The intervention function must return a binary treatment value for every row.")
        as.numeric(values)
    }
    # Evaluate the rule at both possible observed assignments. Equal outputs
    # describe a fixed decision; outputs 0 and 1 preserve the natural value.
    # The remaining case flips treatment and needs different weights.
    value0 <- evaluate(0)
    value1 <- evaluate(1)
    if (any(value0 == 1 & value1 == 0))
        stop("An intervention that flips current treatment requires modified-treatment-policy weights and is not supported.")
    natural <- value0 == 0 & value1 == 1
    list(value = ifelse(natural, observed, value0), natural = natural)
}
