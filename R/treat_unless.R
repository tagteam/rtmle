### treat_unless.R ---
#----------------------------------------------------------------------
## Apply a static intervention unless a contraindication has occurred
#----------------------------------------------------------------------

#' Apply a treatment rule unless a contraindication has occurred
#'
#' This is a convenience intervention function for dynamic treatment rules. It
#' first applies the static values in \code{intervention_table} with
#' \code{\link{intervene}}. For each treatment assignment, it then checks the
#' selected history window of the variables in \code{contra_indication} and
#' changes the assignment to the corresponding value in \code{action} when a
#' contraindication is present.
#'
#' @param data The current history data. Omit this together with
#'   \code{intervention_table} and \code{time_node} to create a configured
#'   intervention function.
#' @param intervention_table Long-format intervention table supplied by
#'   \code{\link{regime}}. Omit this together with \code{data} and
#'   \code{time_node} to create a configured intervention function.
#' @param time_node Current intervention node. Only treatment assignments up
#'   to and including this node are applied. Omit this together with
#'   \code{data} and \code{intervention_table} to create a configured
#'   intervention function.
#' @param contra_indication Character vector of contraindication variable names
#'   without time suffixes, such as \code{c("bleeding", "surgery")}. A row is
#'   contraindicated when any available historical value is non-zero or
#'   \code{TRUE}.
#' @param action Values to assign to the treatment variables after a
#'   contraindication. A single value is recycled across treatment variables;
#'   otherwise values must be supplied in the order of the treatment variables
#'   in \code{intervention_table}. Named values may instead be supplied using
#'   the unsuffixed treatment names; names for treatment variables not yet
#'   assigned at the current node are ignored. Values whose character
#'   representation matches a treatment-option label are converted to a
#'   factor using those options; thus both \code{0} and \code{"0"} match a
#'   treatment option labelled \code{"0"}. \code{NA} always means that the
#'   original observed treatment value is retained, so no intervention is
#'   applied after the contraindication.
#' @param lookback_window Non-negative integer or \code{Inf}. For a treatment
#'   assignment at node \code{k}, contraindications from nodes
#'   \code{max(0, k - lookback_window)} through \code{k} are considered.
#'   Thus, \code{0} considers only the current node and \code{Inf} considers
#'   all available nodes up to and including the current node.
#' @return A copy of \code{data} with the static intervention applied and the
#'   contraindication rule enforced. If \code{data},
#'   \code{intervention_table}, and \code{time_node} are omitted, a configured
#'   intervention function is returned for use as a \code{regime}
#'   \code{intervene_function}.
#' @details Treatment assignments are evaluated separately by intervention
#'   node. At node \code{k}, the contraindication history consists of the
#'   available columns \code{variable_j}, where \code{j} ranges from
#'   \code{max(0, k - lookback_window)} through \code{k}; the current node is
#'   therefore included in the history window. Treatment values before the first
#'   contraindication are left as specified by the static intervention. Calling
#'   \code{treat_unless} with only \code{contra_indication}, \code{action}, and
#'   optionally \code{lookback_window} creates a configured intervention
#'   function, so it can be supplied directly to \code{\link{regime}} without
#'   a wrapper.
#' @seealso \code{\link{intervene}}, \code{\link{regime}}
#' @examples
#' ## Register a configured rule directly as a regime intervention function.
#' ## regime(...,
#' ##     intervene_function = treat_unless(
#' ##         contra_indication = "bleeding",
#' ##         action = 0
#' ##     ))
#' @export
treat_unless <- function(data,
                         intervention_table,
                         time_node,
                         contra_indication,
                         action,
                         lookback_window = Inf) {
    validate_lookback_window <- function(value) {
        if (!is.numeric(value) ||
            length(value) != 1L ||
            is.na(value) ||
            value < 0 ||
            (!is.infinite(value) && value != floor(value))) {
            stop(
                "Argument lookback_window must be a single non-negative ",
                "integer or Inf."
            )
        }
        as.numeric(value)
    }

    if (missing(data) &&
        missing(intervention_table) &&
        missing(time_node)) {
        if (missing(contra_indication) || missing(action)) {
            stop(
                "A configured treat_unless function needs both ",
                "contra_indication and action."
            )
        }
        configured_contra_indication <- force(contra_indication)
        configured_action <- force(action)
        configured_lookback_window <- validate_lookback_window(lookback_window)
        return(function(data, intervention_table, time_node) {
            treat_unless(
                data = data,
                intervention_table = intervention_table,
                time_node = time_node,
                contra_indication = configured_contra_indication,
                action = configured_action,
                lookback_window = configured_lookback_window
            )
        })
    }
    if (missing(data) ||
        missing(intervention_table) ||
        missing(time_node)) {
        stop(
            "data, intervention_table, and time_node must be supplied ",
            "when applying a treat_unless rule."
        )
    }
    if (missing(contra_indication) || missing(action)) {
        stop(
            "contra_indication and action must be supplied when applying ",
            "a treat_unless rule."
        )
    }
    lookback_window <- validate_lookback_window(lookback_window)

    ## Apply the underlying static rule before changing any assignments.
    intervened_data <- data.table::as.data.table(
        intervene(data, intervention_table, time_node)
    )
    observed_data <- data.table::as.data.table(data.table::copy(data))

    if (missing(contra_indication) || is.null(contra_indication)) {
        contra_indication <- character()
    } else if (!is.character(contra_indication) ||
               anyNA(contra_indication) ||
               any(!nzchar(contra_indication)) ||
               any(grepl("_[0-9]+$", contra_indication))) {
        stop(
            "Argument contra_indication must be a character vector of " ,
            "variable names without time suffixes."
        )
    }
    contra_indication <- unique(contra_indication)

    history_table <- data.table::as.data.table(
        data.table::copy(intervention_table)
    )
    if ("time_node" %in% names(history_table)) {
        history_table <- history_table[
            !is.na(history_table[["time_node"]]) &
                history_table[["time_node"]] <= time_node
        ]
    }
    if ("value" %in% names(history_table)) {
        history_table <- history_table[!is.na(history_table[["value"]])]
    }

    if (NROW(history_table) == 0L ||
        !all(c("variable", "value") %in% names(history_table)) ||
        length(contra_indication) == 0L) {
        return(intervened_data)
    }

    treatment_variables <- unique(sub(
        "_[0-9]+$",
        "",
        as.character(history_table[["variable"]])
    ))
    if (length(action) == 0L || is.list(action) ||
        (!is.atomic(action))) {
        stop(
            "Argument action must contain one value or one value for each " ,
            "treatment variable."
        )
    }
    action_names <- names(action)
    if (!is.null(action_names)) {
        if (length(action_names) != length(action) ||
            anyNA(action_names) || any(!nzchar(action_names)) ||
            anyDuplicated(action_names) ||
            any(!treatment_variables %in% action_names)) {
            stop(
                "Named argument action must name each treatment variable " ,
                "without a time suffix."
            )
        }
        action_by_treatment <- action[treatment_variables]
    } else {
        if (length(action) != 1L &&
            length(action) != length(treatment_variables)) {
            stop(
                "Argument action must contain one value or one value for " ,
                "each treatment variable."
            )
        }
        action_by_treatment <- rep(action, length.out = length(treatment_variables))
        names(action_by_treatment) <- treatment_variables
    }

    treatment_options <- attr(
        history_table,
        "treatment_options",
        exact = TRUE
    )
    if (is.null(treatment_options)) {
        treatment_options <- list()
    }
    option_factor <- function(value, treatment_variable) {
        if (length(value) != 1L || is.na(value)) {
            return(value)
        }
        options <- treatment_options[[treatment_variable]]
        if (is.null(options) &&
            is.factor(history_table[["value"]])) {
            options <- levels(history_table[["value"]])
        }
        option_labels <- as.character(options)
        value_label <- as.character(value)
        if (length(option_labels) > 0L &&
            length(value_label) == 1L &&
            !is.na(value_label) &&
            value_label %in% option_labels) {
            return(factor(value_label, levels = option_labels))
        }
        value
    }
    action_names <- names(action_by_treatment)
    action_by_treatment <- lapply(
        seq_along(action_by_treatment),
        function(index) {
            option_factor(
                action_by_treatment[[index]],
                action_names[[index]]
            )
        }
    )
    names(action_by_treatment) <- action_names

    available_names <- names(observed_data)
    has_contra_indication <- vapply(
        contra_indication,
        function(reference) {
            if (reference %in% available_names) return(TRUE)
            suffixed <- available_names[startsWith(
                available_names,
                paste0(reference, "_")
            )]
            any(grepl(
                "^[0-9]+$",
                substring(suffixed, nchar(reference) + 2L)
            ))
        },
        logical(1L)
    )
    if (any(!has_contra_indication)) {
        stop(
            "The following contra_indication variables are not present in " ,
            "data: ",
            paste(contra_indication[!has_contra_indication], collapse = ", ")
        )
    }

    is_contraindicated <- function(values) {
        if (is.logical(values)) {
            return(!is.na(values) & values)
        }
        if (is.numeric(values)) {
            return(!is.na(values) & values != 0)
        }
        as.character(values) %in% c("1", "TRUE", "true")
    }

    coerce_action <- function(value, current) {
        if (is.factor(current)) {
            return(factor(as.character(value), levels = levels(current)))
        }
        if (is.factor(value)) {
            value <- as.character(value)
            if (is.integer(current)) return(as.integer(value))
            if (is.numeric(current)) return(as.numeric(value))
            if (is.logical(current)) {
                return(value %in% c("TRUE", "true", "1"))
            }
        }
        if (is.integer(current) && is.numeric(value)) {
            return(as.integer(value))
        }
        if (is.logical(current) && is.character(value)) {
            return(value %in% c("TRUE", "true", "1"))
        }
        value
    }

    for (row in seq_len(nrow(history_table))) {
        variable <- as.character(history_table[["variable"]][[row]])
        treatment_variable <- sub("_[0-9]+$", "", variable)
        action_value <- action_by_treatment[[treatment_variable]]

        action_node <- if ("time_node" %in% names(history_table)) {
            as.integer(as.character(history_table[["time_node"]][[row]]))
        } else {
            as.integer(time_node)
        }
        history_start <- if (is.infinite(lookback_window)) {
            0L
        } else {
            max(0L, action_node - as.integer(lookback_window))
        }
        history_nodes <- seq.int(history_start, action_node)
        history_variables <- unique(unlist(lapply(
            contra_indication,
            function(reference) {
                c(
                    intersect(reference, available_names),
                    intersect(
                        paste0(reference, "_", history_nodes),
                        available_names
                    )
                )
            }
        ), use.names = FALSE))
        if (length(history_variables) == 0L) next

        contraindicated <- Reduce(
            `|`,
            lapply(
                history_variables,
                function(reference) is_contraindicated(observed_data[[reference]])
            ),
            init = rep(FALSE, nrow(observed_data))
        )
        rows <- which(contraindicated)
        if (length(rows) == 0L) next

        replacement <- if (is.na(action_value)) {
            observed_data[[variable]][rows]
        } else {
            rep(
                coerce_action(action_value, intervened_data[[variable]]),
                length.out = length(rows)
            )
        }
        data.table::set(
            intervened_data,
            i = rows,
            j = variable,
            value = replacement
        )
    }
    intervened_data
}

######################################################################
### treat_unless.R ends here
