#----------------------------------------------------------------------
## Define a treatment regime for an emulated target trial
#----------------------------------------------------------------------

##' Define a treatment regime for an emulated target trial
##'
##' Adds a named regime to an existing \code{rtmle} object. A regime
##' defines treatment values at each time point during follow-up.
##'
##' @param x An \code{rtmle} object as returned by \code{\link{rtmle_init}}.
##' @param name Name of the regime.
##' @param intervention A vector or a data frame (tibble, data table) specifying
##'   factor columns with the treatment values dictated by the regime.
##'   A data frame must contain a time column with values that match all or a
##'   subset of \code{x$time_grid} and columns for all treatments that each contain factors
##'   with two levels allowing \code{NA} to mean no intervention on this variable at that time node.  
##'   The vector version is a shortcut to create the data.frame for regimes that are dictating
##'   the same value for all time nodes. It has 
##'   a single 0/1 value for each element of argument \code{treatment_variables}. 
##' @param expand Logical. If \code{FALSE} and \code{intervention} contains a
##'   time column, do not expand static interventions across time.
##' @param treatment_variables Names of the treatment variable(s) when
##'   \code{intervention} is supplied as a 0/1 vector.
##' @param intervene_function A function, or a character string naming one,
##'   called as \code{fun(data, intervention_table, time_node)}. It must return
##'   the complete intervention-updated data while preserving its rows and
##'   required columns. For a dynamic treatment regime, use it to update
##'   treatment values according to the observed history.
##' @param intervene_function_args Named list of additional arguments passed to
##'   \code{intervene_function} at every intervention node. This is useful for
##'   parameterized functions such as \code{\link{treat_unless}}, for example
##'   \code{list(contra_indication = "bleeding", action = 0)}. The names
##'   \code{data}, \code{intervention_table}, and the time argument are reserved.
##' @param multiple_treatment_factorization How propensity models for multiple
##'   treatment variables assigned at the same node are factorized. One of
##'   \code{"joint"}, \code{"sequential"}, or \code{"independent"}; the
##'   default is \code{"joint"}.
##' @param adherence_model_strata Optional character vector of prepared-data
##'   variables defining subgroups in which separate adherence propensities are
##'   fitted. This argument is only for dynamic treatment rules. If it is
##'   \code{NULL}, one pooled binary adherence model is fitted. Unsuffixed
##'   time-varying names refer to the latest pre-decision value. At a discrete
##'   node \code{k}, \code{bleeding} resolves to \code{bleeding_{k-1}}
##'   (or \code{bleeding_0} at baseline). At a continuous treatment event,
##'   it resolves to \code{bleeding_k}, including measurements made at the
##'   treatment event. The current treatment assignment is not pre-decision.
##'   The argument defines fitting strata and does not add model predictors.
##'   Multiple variables define strata by their combinations.
##' @param verbose Logical. If \code{FALSE}, suppress all messages.
##' @param ... Additional arguments are not used. The former
##'   \code{dynamic_propensity_instructions} interface is no longer supported;
##'   use \code{adherence_model_strata} instead.
##' @return The modified object contains the treatment variables and
##'   \code{intervention_table} in \code{x$regimes[[name]]}.
##' @details
##'   For a dynamic rule, \code{rtmle} constructs a binary adherence response
##'   indicating whether the observed treatment agrees with the
##'   intervention-updated treatment. \code{adherence_model_strata} controls
##'   whether this response is modelled in one pooled model or in separate
##'   pre-decision subgroups. The intervention function determines the regime
##'   itself, including whether a contraindication pauses treatment for one
##'   interval or stops it permanently. A static regime uses the ordinary
##'   observed-treatment model. Additional intervention-function arguments can
##'   be supplied through \code{intervene_function_args}, avoiding a wrapper
##'   when a dynamic rule such as \code{\link{treat_unless}} has fixed
##'   parameters.
##'
##'   In continuous time, supply a single binary nominal value for treatment
##'   \code{A} and optionally an intervention function. The function is
##'   evaluated at baseline and at observed treatment opportunities; its
##'   \code{time_node} is the event index, not calendar time. It receives the
##'   history through that event and an intervention table containing the
##'   nominal value at each index up to that event. It must return a binary
##'   \code{A_<time_node>} for each subject. Covariate-only events contribute
##'   to the history but do not themselves assign treatment. The optional
##'   \code{adherence_model_strata} specifies separate adherence models using
##'   variables known before each treatment decision. A rule may retain the observed treatment
##'   for some histories, with adherence probability one for those histories.
##'   Rules that flip the current observed treatment are not supported, as
##'   they require different weighting.
##' @seealso \code{\link{rtmle_init}}, \code{\link{prepare_rtmle_data}},
##'   \code{\link{intervention_match}}, \code{\link{target}},
##'   \code{\link{model_formula}}, \code{\link{run_rtmle}}
##' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
##' @examples
##' x <- rtmle_init(time_grid = 0:3, name_id = "id",
##'                 name_outcome = "Y", name_competing = "Dead",
##'                 name_censoring = "Censored", censored_label = "censored")
##' x <- regime(
##'     x,
##'     name = "Always_A_never_B",
##'     intervention = data.frame(
##'         time_node = x$intervention_nodes,
##'         A = factor("1", levels = c("0", "1")),
##'         B = factor("0", levels = c("0", "1"))
##'     ),
##'     multiple_treatment_factorization = "sequential"
##' )
##'
##' # A dynamic rule can stop A after a contraindication in the history.
##' x <- regime(
##'     x, name = "A_until_bleeding",
##'     intervention = data.frame(
##'         time_node = x$intervention_nodes,
##'         A = factor("1", levels = c("0", "1"))
##'     ),
##'     intervene_function = treat_unless(
##'         contra_indication = "bleeding",
##'         action = 0
##'     ),
##'     adherence_model_strata = "bleeding",
##'     verbose = FALSE
##' )
##' @export
protocol <- function(x,
                     name,
                     intervention,
                     expand = TRUE,
                     treatment_variables,
                     intervene_function = NULL,
                     verbose = TRUE,
                     multiple_treatment_factorization = "joint",
                     adherence_model_strata = NULL,
                     intervene_function_args = NULL,
                     ...) {
    if (isTRUE(x$continuous_time)) {
        if (length(list(...))) stop("Unused arguments for continuous-time regime.")
        validate_adherence_model_strata(adherence_model_strata)
        if (length(adherence_model_strata) && is.null(intervene_function))
            stop("`adherence_model_strata` is only for dynamic treatment rules; supply `intervene_function`.")
        if (!identical(multiple_treatment_factorization, "joint"))
            stop("Continuous time supports treatment A only.")
        if (!is.null(intervene_function) &&
            !(is.function(intervene_function) ||
              (is.character(intervene_function) && length(intervene_function) == 1L &&
               !is.na(intervene_function) && nzchar(intervene_function))))
            stop("intervene_function must be a function or the name of a function.")
        if (!is.null(intervene_function_args) && !is.list(intervene_function_args))
            stop("intervene_function_args must be a named list.")
        if (length(intervene_function_args)) {
            nms <- names(intervene_function_args)
            if (is.null(intervene_function) || is.null(nms) || anyNA(nms) ||
                any(!nzchar(nms)) || anyDuplicated(nms) ||
                any(nms %in% c("data", "intervention_table", "time_node", "time", "current_time", "current.time")))
                stop("intervene_function_args requires an intervention function and unique, non-reserved argument names.")
        }
        x <- continuous_regime(x, name, intervention,
                               if (missing(treatment_variables)) "A" else treatment_variables)
        if (!is.null(intervene_function)) {
            x$regimes[[name]]$dynamic_intervention <- TRUE
            x$regimes[[name]]$intervene_function <- intervene_function
            x$regimes[[name]]$intervene_function_args <- intervene_function_args
            x$regimes[[name]]$adherence_model_strata <- adherence_model_strata
        }
        return(x)
    }
    dots <- list(...)
    if (length(dots) > 0L) {
        dot_names <- names(dots)
        if (!is.null(dot_names) &&
            "dynamic_propensity_instructions" %in% dot_names) {
            stop(
                "The former interface was removed: ",
                "dynamic_propensity_instructions; use adherence_model_strata instead."
            )
        }
        stop(
            "Unused argument(s): ",
            paste(dot_names[!is.na(dot_names) & nzchar(dot_names)], collapse = ", ")
        )
    }
    variable <- time_node <- NULL
    if (length(multiple_treatment_factorization) != 1L || !is.character(multiple_treatment_factorization) || is.na(multiple_treatment_factorization) ||
        !multiple_treatment_factorization %in% c("joint", "sequential", "independent")) {
        stop("Argument multiple_treatment_factorization must be one of: ", paste(c("joint", "sequential", "independent"), collapse = ", "), ".")
    }
    validate_adherence_model_strata(adherence_model_strata, "Argument adherence_model_strata")
    if (length(adherence_model_strata) == 0L) {
        adherence_model_strata <- NULL
    }
    if (!is.null(adherence_model_strata) && length(intervene_function) == 0L) {
        stop(
            "`adherence_model_strata` is only for dynamic treatment rules; ",
            "supply `intervene_function`."
        )
    }

    allowed_intervention_node_names <- c(
        "time_node", "intervention_node", "node", "time_grid", "time"
    )
    if (inherits(intervention, "data.frame")) {
        intervention_table <- data.table::copy(intervention)
        data.table::setDT(intervention_table)
        treatment_variables <- names(intervention_table)
        intervention_node_name <- intersect(
            treatment_variables,
            allowed_intervention_node_names
        )
        if (length(intervention_node_name) > 0L) {
            if (length(intervention_node_name) > 1L) {
                intervention_node_name <- intervention_node_name[[1L]]
            }
            treatment_variables <- setdiff(
                treatment_variables,
                intervention_node_name
            )
        } else if (length(x$intervention_nodes) > 1L) {
            stop(
                "Argument intervention needs to have a variable called ",
                "'intervention_node' with values equal to or a subset of ",
                "x$intervention_nodes."
            )
        }
        if (any(grepl("_[0-9]+$", treatment_variables))) {
            stop("Treatment variables should be given without time suffix.")
        }
        if (any(is.na(match(
            intervention_table[["time_node"]],
            x$intervention_nodes,
            nomatch = NA
        )))) {
            too_many <- which(is.na(match(
                intervention_table[["time_node"]],
                x$intervention_nodes,
                nomatch = NA
            )))
            stop(
                "The following time points are not registered as intervention ",
                "nodes in the object: ",
                paste(intervention_table[["time_node"]][too_many], collapse = ", ")
            )
        }
        data.table::setnames(
            intervention_table,
            old = intervention_node_name,
            new = "time_node"
        )
        treatment_options <- lapply(treatment_variables, function(variable) {
            if (!is.factor(intervention[[variable]])) {
                stop(
                    "The treatment variables must be factors. Problem with ",
                    variable, "."
                )
            }
            if (length(levels(intervention[[variable]])) != 2L) {
                stop(
                    "All treatment variables must have exactly 2 levels. ",
                    "Problem with variable ", variable, "."
                )
            }
            levels(intervention[[variable]])
        })
        names(treatment_options) <- treatment_variables
    } else {
        if (!missing(treatment_variables) &&
            length(treatment_variables) == length(intervention) &&
            all(intervention %in% c(0, 1))) {
            intervention_table <- data.table::as.data.table(lapply(
                seq_along(treatment_variables),
                function(index) factor(intervention[[index]], levels = c(0, 1))
            ))
            intervention_table <- cbind(
                time_node = x$intervention_nodes,
                intervention_table
            )
            data.table::setnames(
                intervention_table,
                c("time_node", treatment_variables)
            )
            treatment_options <- lapply(treatment_variables, function(x) c(0, 1))
            names(treatment_options) <- treatment_variables
        } else {
            stop(
                "Argument intervention is not a data.frame. Hence it must ",
                "be a vector of 0s and 1s with the same length as ",
                "argument treatment_variables."
            )
        }
    }

    intervention_table <- data.table::melt(
        intervention_table,
        id.vars = "time_node",
        variable.name = "variable",
        value.name = "value",
        value.factor = TRUE
    )[, variable := paste0(variable, "_", time_node)]

    if (length(intervene_function) > 0L &&
        !(is.function(intervene_function) ||
          (is.character(intervene_function) && length(intervene_function) == 1L &&
           !is.na(intervene_function) && nzchar(intervene_function)))) {
        stop("intervene_function must be a function or the name of a function.")
    }
    if (!is.null(intervene_function_args) &&
        !is.list(intervene_function_args)) {
        stop("intervene_function_args must be a named list.")
    }
    if (length(intervene_function_args) > 0L) {
        argument_names <- names(intervene_function_args)
        reserved_names <- c(
            "data",
            "intervention_table",
            "time_node",
            "time",
            "current_time",
            "current.time"
        )
        if (is.null(argument_names) ||
            anyNA(argument_names) ||
            any(!nzchar(argument_names)) ||
            anyDuplicated(argument_names) ||
            any(argument_names %in% reserved_names)) {
            stop(
                "intervene_function_args must be a named list with unique ",
                "non-reserved argument names."
            )
        }
        if (length(intervene_function) == 0L) {
            stop(
                "intervene_function_args requires an intervene_function."
            )
        }
    }

    x$regimes[[name]] <- NULL
    x$regimes[[name]]$dynamic_intervention <- length(intervene_function) > 0L
    x$regimes[[name]]$intervene_function <- if (length(intervene_function) > 0L) {
        intervene_function
    } else {
        "intervene"
    }
    x$regimes[[name]]$intervene_function_args <- intervene_function_args
    x$regimes[[name]]$treatment_variables <- treatment_variables
    x$regimes[[name]]$treatment_options <- treatment_options
    x$regimes[[name]]$intervention_table <- intervention_table[]
    x$regimes[[name]]$multiple_treatment_factorization <-
        multiple_treatment_factorization
    x$regimes[[name]]$adherence_model_strata <- adherence_model_strata

    if (length(x$names$treatment_options) == 0L) {
        x$names$treatment_options <- treatment_options
    } else {
        new_options <- setdiff(
            names(treatment_options),
            names(x$names$treatment_options)
        )
        if (length(new_options) > 0L) {
            x$names$treatment_options <- c(
                x$names$treatment_options,
                treatment_options[new_options]
            )
        }
    }
    x <- intervention_match(x, regime_name = name)
    x
}

######################################################################
##' @rdname protocol
##' @export
regime <- protocol

######################################################################
### protocol.R ends here

validate_adherence_model_strata <- function(strata, label = "`adherence_model_strata`") {
    if (!is.null(strata) &&
        (!is.character(strata) || anyNA(strata) || any(!nzchar(strata)))) {
        stop(label, " must be NULL or a character vector of variable names.")
    }
    invisible(NULL)
}
