### evaluate_intervention.R ---
#----------------------------------------------------------------------
## Evaluate and validate a user-defined intervention function
#----------------------------------------------------------------------

evaluate_intervention <- function(regime,
                                  data,
                                  intervention_table,
                                  time_node,
                                  full_n = NROW(data),
                                  row_indices = seq_len(NROW(data))) {
    if (length(full_n) != 1L || is.na(full_n) || !is.finite(full_n) ||
        full_n < 0L || full_n != as.integer(full_n)) {
        stop("`full_n` must be one non-negative integer.")
    }
    if (!is.numeric(row_indices) || length(row_indices) != NROW(data) ||
        anyNA(row_indices) || any(row_indices < 1L) ||
        any(row_indices > full_n)) {
        stop("`row_indices` must identify the rows of `data` in the full cohort.")
    }
    if (length(regime$intervene_function) == 0) {
        stop("The regime has no intervene_function.")
    }

    history_table <- data.table::copy(intervention_table)
    if (NROW(history_table) > 0L && "time_node" %in% names(history_table)) {
        keep <- !is.na(history_table[["time_node"]]) &
            history_table[["time_node"]] <= time_node
        if ("value" %in% names(history_table)) {
            keep <- keep & !is.na(history_table[["value"]])
        }
        history_table <- history_table[keep]
    }
    if (!is.null(regime$treatment_options)) {
        data.table::setattr(
            history_table,
            "treatment_options",
            regime$treatment_options
        )
    }

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
        function_args <- regime$intervene_function_args
        if (is.null(function_args)) {
            function_args <- list()
        }
        if (!is.list(function_args) ||
            (length(function_args) > 0L &&
             (is.null(names(function_args)) ||
              anyNA(names(function_args)) ||
              any(!nzchar(names(function_args))) ||
              anyDuplicated(names(function_args))))) {
            stop("intervene_function_args must be a named list.")
        }
        args <- c(list(
            data = data.table::copy(data),
            intervention_table = data.table::copy(history_table)
        ), function_args)
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

    intervened_data <- call_with_context(regime$intervene_function)
    if (!(inherits(intervened_data, "data.frame") ||
          is.matrix(intervened_data))) {
        stop(
            "The intervene_function must return the intervention-updated data ",
            "as a data frame or matrix; do not return adherence-model metadata."
        )
    }
    if (NROW(intervened_data) != NROW(data)) {
        stop("The intervention function must preserve the number and order of rows.")
    }
    missing_columns <- setdiff(names(data), colnames(intervened_data))
    if (length(missing_columns) > 0L) {
        stop(
            "The intervention function removed required data column(s): ",
            paste(missing_columns, collapse = ", "), "."
        )
    }
    intervened_data <- data.table::as.data.table(intervened_data)

    validate_adherence_model_strata(regime$adherence_model_strata)
    adherence_model_strata <- resolve_adherence_model_strata(
        regime$adherence_model_strata,
        time_node = time_node,
        available_names = names(data)
    )
    missing_strata <- setdiff(adherence_model_strata, names(data))
    if (length(missing_strata) > 0L) {
        stop(
            "Unknown variable(s) supplied in `adherence_model_strata`: ",
            paste(missing_strata, collapse = ", ")
        )
    }

    list(
        data = intervened_data,
        dynamic_adherence = isTRUE(regime$dynamic_intervention),
        adherence_model_strata = adherence_model_strata,
        multiple_treatment_factorization =
            regime$multiple_treatment_factorization
    )
}

######################################################################
### evaluate_intervention.R ends here
