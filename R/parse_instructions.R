### parse_instructions.R ---
#----------------------------------------------------------------------
## Validate adherence-model stratification supplied through regime()
#----------------------------------------------------------------------

validate_adherence_model_strata <- function(
    strata,
    label = "`adherence_model_strata`"
) {
    if (is.null(strata)) {
        return(invisible(NULL))
    }
    if (!is.character(strata) || anyNA(strata) || any(!nzchar(strata))) {
        stop(label, " must be NULL or a character vector of variable names.")
    }
    invisible(NULL)
}

resolve_adherence_model_strata <- function(strata, time_node, available_names) {
    validate_adherence_model_strata(strata)
    if (is.null(strata) || length(strata) == 0L) {
        return(NULL)
    }
    if (length(time_node) != 1L || is.na(time_node) ||
        !is.finite(time_node)) {
        stop("`time_node` must be one finite value.")
    }
    time_node <- as.integer(time_node)
    resolved <- vapply(strata, function(reference) {
        # Explicit time suffixes are useful for lagged variables, but a
        # current- or future-node variable would violate the pre-decision
        # history restriction. Node zero is the baseline/pre-first-decision
        # convention used throughout the prepared data.
        if (grepl("_[0-9]+$", reference)) {
            suffix <- as.integer(sub(".*_", "", reference))
            if (suffix > 0L && suffix >= time_node) {
                stop(
                    "`adherence_model_strata` must use variables known before ",
                    "the treatment decision at node ", time_node, ". `",
                    reference, "` is not pre-decision there."
                )
            }
            return(reference)
        }

        # An unsuffixed time-varying name means the latest available value,
        # not the value from the current interval.
        lagged_reference <- paste0(reference, "_", max(0L, time_node - 1L))
        if (lagged_reference %in% available_names) {
            return(lagged_reference)
        }
        reference
    }, character(1))
    unique(resolved)
}

validate_multiple_treatment_factorization <- function(value, label) {
    allowed <- c("joint", "sequential", "independent")
    if (length(value) != 1L || !is.character(value) || is.na(value) ||
        !value %in% allowed) {
        stop(label, " must be one of: ", paste(allowed, collapse = ", "), ".")
    }
    invisible(value)
}

######################################################################
### parse_instructions.R ends here
