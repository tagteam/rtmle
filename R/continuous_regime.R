continuous_regime <- function(x, name, intervention, treatment_variables) {
    if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) stop("Regime name must be nonempty.")
    if (is.data.frame(intervention)) {
        if (!"A" %in% names(intervention)) stop("Continuous time supports treatment A only.")
        if (length(setdiff(names(intervention), c("A", "time", "time_node", "node", "time_grid", "intervention_node"))))
            stop("Continuous time supports treatment A only.")
        time_columns <- setdiff(names(intervention), "A")
        if (nrow(intervention) != 1L || any(vapply(intervention[time_columns], function(t) {
            anyNA(t) || !all(as.character(t) == "0")
        }, logical(1))))
            stop("Continuous time requires a single static treatment value at baseline; scheduled interventions are not supported.")
        intervention <- unique(as.character(intervention$A))
    }
    if (!identical(treatment_variables, "A") || length(intervention) != 1L ||
        is.na(intervention) || !intervention %in% c(0, 1))
        stop("Continuous time supports a static intervention of 0 or 1 on A only.")
    x$regimes[[name]] <- list(static_intervention = as.numeric(intervention),
                              treatment_variables = "A", treatment_options = list(A = c("0", "1")),
                              dynamic_intervention = FALSE,
                              intervention_table = data.table::data.table(time_node = 0, variable = "A_0", value = as.character(intervention)))
    x$names$treatment_options <- list(A = c("0", "1"))
    x$estimate <- x$IC <- NULL
    x$continuous_nuisance[[name]] <- NULL
    x$sequential_outcome_regression <- x$run_time_horizons <- NULL
    x
}
