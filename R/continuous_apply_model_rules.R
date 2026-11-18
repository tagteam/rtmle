# Apply the same named exclusion and inclusion rules used by model_formula.
# `available` contains terms observed by the time the model is fitted.
continuous_apply_model_rules <- function(history, available, all_names, model_options,
                                         outcome_name, time_covariates) {
    excluded <- model_options$exclude_variables
    # Custom marginal censoring formulas can contain terms such as
    # strata(group) or interactions. Excluding a variable removes any such
    # term that depends on it as well.
    if (length(excluded) && length(history)) {
        is_excluded <- function(variables) {
            matches <- variables %in% excluded
            for (variable in intersect(excluded, time_covariates)) {
                prefix <- paste0(variable, "_")
                suffix <- substring(variables, nchar(prefix) + 1L)
                matches <- matches | (startsWith(variables, prefix) & grepl("^[0-9]+$", suffix))
            }
            any(matches)
        }
        drop <- vapply(history, function(term) {
            variables <- tryCatch(all.vars(stats::as.formula(paste("~", term))),
                                  error = function(e) character())
            is_excluded(c(term, variables))
        }, logical(1))
        history <- history[!drop]
    }

    matching_rules <- function(rules) {
        if (!length(rules) || is.null(outcome_name)) return(character())
        matches <- vapply(names(rules), grepl, logical(1), x = outcome_name)
        unlist(rules[matches], use.names = FALSE)
    }
    patterns <- matching_rules(model_options$exclusion_rules)
    if (length(patterns) && length(history)) {
        remove <- vapply(history, function(name) any(vapply(patterns, grepl, logical(1), x = name)), logical(1))
        history <- history[!remove]
    }
    for (pattern in matching_rules(model_options$inclusion_rules)) {
        additions <- grep(pattern, available, value = TRUE)
        if (!length(additions)) {
            if (length(grep(pattern, all_names, value = TRUE)))
                stop("Inclusion rule for ", outcome_name,
                     " matches only variables unavailable at that event: ", pattern, ".")
            warning("Inclusion rule for ", outcome_name,
                    " did not match an available variable: ", pattern, ".", call. = FALSE)
        }
        history <- unique(c(history, additions))
    }
    history
}
