### parse_instructions.R ---
#----------------------------------------------------------------------
## Parse propensity instructions supplied through protocol()
#----------------------------------------------------------------------

parse_instructions <- function(propensity_instructions,
                               n,
                               full_n = n,
                               allow_full_length = FALSE) {
    parse_propensity_instruction <- function(spec, label) {
        if (is.numeric(spec) && is.null(dim(spec))) {
            probability <- spec
            if (length(probability) == 1L && n != 1L) {
                probability <- rep(probability, n)
            }
            if (length(probability) != n &&
                !(allow_full_length && length(probability) == full_n)) {
                stop("The `probability` in ", label,
                     " must have length 1 or NROW(data).")
            }
            if (any(is.nan(probability))) {
                stop("The `probability` in ", label,
                     " cannot contain NaN; use NA to request estimation.")
            }
            observed <- probability[!is.na(probability)]
            if (any(!is.finite(observed)) ||
                any(observed < 0 | observed > 1)) {
                stop("Non-missing probabilities in ", label,
                     " must be finite values in [0, 1].")
            }
            return(list(
                mode = "fixed",
                probability = probability,
                stratify_by = NULL
            ))
        }
        if (!is.list(spec) || is.null(spec$mode)) {
            stop(label, " must be a numeric fixed-probability vector or a list with `mode`.")
        }
        mode <- as.character(spec$mode)
        if (length(mode) != 1L || is.na(mode) ||
            !mode %in% c("fixed", "adherence")) {
            stop("The `mode` in ", label,
                 " must be `fixed` or `adherence`.")
        }
        unknown <- setdiff(names(spec), c("mode", "probability", "stratify_by"))
        if (length(unknown) > 0) {
            stop("Unknown element(s) in ", label, ": ",
                 paste(unknown, collapse = ", "), ".")
        }
        stratify_by <- spec$stratify_by
        if (length(stratify_by) > 0) {
            if (!is.character(stratify_by) || anyNA(stratify_by) ||
                any(!nzchar(stratify_by))) {
                stop("`stratify_by` in ", label,
                     " must contain non-missing variable names.")
            }
            stratify_by <- unique(stratify_by)
        } else {
            stratify_by <- NULL
        }
        if (identical(mode, "fixed")) {
            if (length(stratify_by) > 0) {
                stop("`stratify_by` can only be used with `mode = \"adherence\"`.")
            }
            if (is.null(spec$probability)) {
                stop("A fixed propensity instruction must supply `probability`.")
            }
            probability <- spec$probability
            if (!is.numeric(probability) || !is.null(dim(probability))) {
                stop("The `probability` in ", label,
                     " must be a numeric vector.")
            }
            if (length(probability) == 1L && n != 1L) {
                probability <- rep(probability, n)
            }
            if (length(probability) != n &&
                !(allow_full_length && length(probability) == full_n)) {
                stop("The `probability` in ", label,
                     " must have length 1 or NROW(data).")
            }
            if (any(is.nan(probability))) {
                stop("The `probability` in ", label,
                     " cannot contain NaN; use NA to request estimation.")
            }
            observed <- probability[!is.na(probability)]
            if (any(!is.finite(observed)) ||
                any(observed < 0 | observed > 1)) {
                stop("Non-missing probabilities in ", label,
                     " must be finite values in [0, 1].")
            }
        } else {
            if (!is.null(spec$probability)) {
                stop("An adherence propensity instruction must not supply `probability`.")
            }
            probability <- NULL
        }
        list(
            mode = mode,
            probability = probability,
            stratify_by = stratify_by
        )
    }

    if (length(propensity_instructions) == 0L) {
        return(NULL)
    }
    if (is.numeric(propensity_instructions) &&
        is.null(dim(propensity_instructions))) {
        return(list(
            .default = parse_propensity_instruction(
                propensity_instructions,
                "the unnamed propensity instruction"
            )
        ))
    }
    if (inherits(propensity_instructions, "data.frame") ||
        is.matrix(propensity_instructions)) {
        instruction_names <- colnames(propensity_instructions)
        if (is.null(instruction_names) || any(!nzchar(instruction_names))) {
            stop("A table of propensity instructions must have treatment-variable column names.")
        }
        parsed <- lapply(instruction_names, function(variable) {
            parse_propensity_instruction(
                propensity_instructions[[variable]],
                paste0("propensity instruction `", variable, "`")
            )
        })
        names(parsed) <- instruction_names
        return(parsed)
    }
    if (is.list(propensity_instructions)) {
        instruction_names <- names(propensity_instructions)
        if (!is.null(propensity_instructions$mode)) {
            return(list(
                .default = parse_propensity_instruction(
                    propensity_instructions,
                    "the default propensity instruction"
                )
            ))
        }
        if (is.null(instruction_names) || any(!nzchar(instruction_names))) {
            stop("Propensity instructions must be named by treatment variable.")
        }
        parsed <- lapply(instruction_names, function(variable) {
            parse_propensity_instruction(
                propensity_instructions[[variable]],
                paste0("propensity instruction `", variable, "`")
            )
        })
        names(parsed) <- instruction_names
        return(parsed)
    }
    stop("Unsupported `propensity_instructions`; use fixed numeric vectors or named instruction lists.")
}

######################################################################
### parse_instructions.R ends here
