### evaluate_intervention.R ---
#----------------------------------------------------------------------
## Evaluate and validate a user-defined intervention function
#----------------------------------------------------------------------

evaluate_intervention <- function(protocol,
                                  data,
                                  intervention_table,
                                  time_node,
                                  full_n = NROW(data),
                                  row_indices = seq_len(NROW(data))) {
    if (length(full_n) != 1L || is.na(full_n) || !is.finite(full_n) ||
        full_n < 0L ||
        full_n != as.integer(full_n)) {
        stop("`full_n` must be one non-negative integer.")
    }
    if (!is.numeric(row_indices) || length(row_indices) != NROW(data) ||
        anyNA(row_indices) || any(row_indices < 1L) ||
        any(row_indices > full_n)) {
        stop("`row_indices` must identify the rows of `data` in the full cohort.")
    }
    if (length(protocol$intervene_function) == 0) {
        stop("The protocol has no intervene_function.")
    }
    history_table <- data.table::copy(intervention_table)
    if (NROW(history_table) > 0 && "time_node" %in% names(history_table)) {
        keep <- !is.na(history_table[["time_node"]]) &
            history_table[["time_node"]] <= time_node
        if ("value" %in% names(history_table)) {
            keep <- keep & !is.na(history_table[["value"]])
        }
        history_table <- history_table[keep]
    }

    # The intervention and optional propensity callbacks share the same
    # stage-specific data, intervention history, and node arguments.
    call_with_context <- function(fun) {
        if (!is.function(fun)) {
            if (!is.character(fun) || length(fun) != 1L || is.na(fun) ||
                !nzchar(fun)) {
                stop("intervene_function must be a function or the name of a function.")
            }
            if (exists(fun, envir = .GlobalEnv, mode = "function", inherits = TRUE)) {
                fun <- get(fun, envir = .GlobalEnv, mode = "function", inherits = TRUE)
            } else {
                fun <- match.fun(fun)
            }
        }
        formal_names <- names(formals(fun))
        args <- list(
            data = data.table::copy(data),
            intervention_table = data.table::copy(history_table)
        )
        recognized_time_names <- c(
            "time_node", "time", "current_time", "current.time"
        )
        time_name <- intersect(recognized_time_names, formal_names)
        if (length(time_name) > 0L) {
            args[[time_name[[1L]]]] <- time_node
        } else if ("..." %in% formal_names) {
            args$time_node <- time_node
        } else if (length(formal_names) >= 3L) {
            args[[formal_names[[3L]]]] <- time_node
        }
        do.call(fun, args)
    }

    intervened_data <- call_with_context(protocol$intervene_function)
    if (!(inherits(intervened_data, "data.frame") ||
          is.matrix(intervened_data))) {
        stop(
            "The intervene_function must return the intervention-updated data ",
            "as a data frame or matrix. Supply propensity_instructions and ",
            "propensity_variables to protocol() instead of returning metadata."
        )
    }

    instructions_are_static <- !is.function(protocol$propensity_instructions)
    propensity_instructions <- protocol$propensity_instructions
    if (is.function(propensity_instructions)) {
        propensity_instructions <- call_with_context(propensity_instructions)
    }
    propensity_instructions <- parse_instructions(
        propensity_instructions,
        NROW(data),
        full_n = full_n,
        allow_full_length = instructions_are_static
    )

    # Static fixed probabilities are often declared once for the full cohort,
    # while nuisance fits use an at-risk subset. Align those vectors with the
    # rows being evaluated before handing them to downstream code. A
    # function-valued instruction is evaluated separately on each subset and
    # therefore must already have the current length.
    if (instructions_are_static && length(propensity_instructions) > 0L &&
        full_n != NROW(data)) {
        for (instruction_name in names(propensity_instructions)) {
            instruction <- propensity_instructions[[instruction_name]]
            if (identical(instruction$mode, "fixed") &&
                length(instruction$probability) == full_n) {
                instruction$probability <- instruction$probability[row_indices]
                propensity_instructions[[instruction_name]] <- instruction
            }
        }
    }

    propensity_variables <- protocol$propensity_variables
    if (is.function(propensity_variables)) {
        propensity_variables <- call_with_context(propensity_variables)
    }
    if (length(propensity_variables) > 0L) {
        if (is.character(propensity_variables)) {
            if (anyNA(propensity_variables) || any(!nzchar(propensity_variables))) {
                stop("`propensity_variables` must contain non-missing prepared-data column names.")
            }
            propensity_variables <- unique(propensity_variables)
        } else if (is.list(propensity_variables) &&
                   !is.null(names(propensity_variables))) {
            valid <- vapply(propensity_variables, function(value) {
                is.character(value) && !anyNA(value) && all(nzchar(value))
            }, logical(1))
            if (!all(valid)) {
                stop("Each `propensity_variables` list element must be a character vector.")
            }
            propensity_variables <- lapply(propensity_variables, unique)
        } else {
            stop("`propensity_variables` must be a character vector or named list.")
        }
    } else {
        propensity_variables <- NULL
    }

    if (NROW(intervened_data) != NROW(data)) {
        stop("The intervention function must preserve the number and order of rows.")
    }
    missing_columns <- setdiff(names(data), colnames(intervened_data))
    if (length(missing_columns) > 0L) {
        stop("The intervention function removed required data column(s): ",
             paste(missing_columns, collapse = ", "), ".")
    }
    intervened_data <- data.table::as.data.table(intervened_data)

    declared_variables <- if (is.list(propensity_variables)) {
        unique(unlist(propensity_variables, use.names = FALSE))
    } else {
        propensity_variables
    }
    missing_propensity_variables <- setdiff(declared_variables, names(data))
    if (length(missing_propensity_variables) > 0L) {
        stop("Unknown propensity_variables supplied to protocol(): ",
             paste(missing_propensity_variables, collapse = ", "), ".")
    }

    if (length(propensity_instructions) > 0L) {
        current_variables <- history_table[["variable"]][
            history_table[["time_node"]] == time_node
        ]
        # Static metadata may declare one instruction for every treatment
        # node, even though the callback is being evaluated at an earlier
        # node. Validate names against the complete protocol table; downstream
        # task selection still uses only the current treatment columns.
        all_intervention_variables <- if (!is.null(protocol$intervention_table) &&
                                          "variable" %in% names(protocol$intervention_table)) {
            protocol$intervention_table[["variable"]]
        } else {
            intervention_table[["variable"]]
        }
        permitted_names <- unique(c(
            all_intervention_variables,
            paste(current_variables, collapse = ","),
            ".default"
        ))
        unknown_instructions <- setdiff(
            names(propensity_instructions),
            permitted_names
        )
        if (length(unknown_instructions) > 0L) {
            stop("Unknown propensity-instruction treatment column(s): ",
                 paste(unknown_instructions, collapse = ", "), ".")
        }
        declared_strata <- unique(unlist(lapply(
            propensity_instructions,
            function(instruction) instruction$stratify_by
        ), use.names = FALSE))
        missing_strata <- setdiff(declared_strata, names(data))
        if (length(missing_strata) > 0L) {
            stop("Unknown `stratify_by` variable(s) supplied to protocol(): ",
                 paste(missing_strata, collapse = ", "), ".")
        }
    }

    list(
        data = intervened_data,
        propensity_instructions = propensity_instructions,
        propensity_variables = propensity_variables
    )
}

######################################################################
### evaluate_intervention.R ends here
