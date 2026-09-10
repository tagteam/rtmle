### model_formula.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 16 2025 (08:58) 
## Version: 
## Last-Updated: maj 21 2026 (08:32) 
##           By: Thomas Alexander Gerds
##     Update #: 162
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Specify formulas for nuisance-parameter models
##' 
##'
##' Because the data are discretized, time-dependent covariates at time \code{k}
##' are not included in treatment models at time \code{k}, except at time 0.
##' This avoids encoding an ordering assumption in which \code{L_k} precedes
##' \code{A_k} when, in the original event history, \code{A_k} may occur before
##' \code{L_k}. A history-dependent intervention can explicitly declare
##' \code{L_k} in the \code{propensity_variables} supplied to
##' \code{\link{protocol}}. Doing so asserts that \code{L_k} is known before
##' the treatment decision at node \code{k}.
##'
##' @title Model formulas for nuisance parameters
##' @param x An object of class \code{"rtmle"} with prepared data and at least
##'   one protocol defined by \code{\link{protocol}}.
##' @param propensity_model Relevant only for interventions depending on
##'   multiple treatment variables. Controls how to model the propensity of the
##'   regimen dictated by the protocol. Possible values are \code{"joint"},
##'   \code{"sequential"}, and \code{"independent"}.
##' @param Markov Names of time-dependent variables that should appear only with
##'   their most recent values on the right-hand side of formulas.
##' @param exclude_variables Variables to exclude from the formulas for the nuisance parameters.
##' @param exclusion_rules Experimental. Additional exclusion rules given as a
##'   named list. Names are variables occurring on the left-hand side of a
##'   formula, and elements are variables to exclude from the right-hand side.
##' @param inclusion_rules Experimental. Additional inclusion rules given as a
##'   named list. Names are variables occurring on the left-hand side of a
##'   formula, and elements are variables to include on the right-hand side.
##' @param verbose Logical. If \code{FALSE} suppress all messages. \code{TRUE} is the default.
##' @param ... Not used.
##' @return The modified \code{rtmle} object.
##' @seealso \code{\link{prepare_rtmle_data}}, \code{\link{protocol}},
##'   \code{\link{target}}, \code{\link{make_exclusion_rule}},
##'   \code{\link{run_rtmle}}
##' @examples
##' data(simulated_cohort)
##' ld <- register_format(simulated_cohort)
##' x <- rtmle_init(time_grid = seq(0,20,4),name_id = "id",name_outcome = "stroke",
##'                 name_competing = "death",
##'                 name_censoring = "dropout",censored_label = "censored")
##' x <- add_long_data(x,
##'                    outcome_data=ld$timevar_data$stroke[!duplicated(id)],
##'                    censored_data=ld$timevar_data$dropout,
##'                    competing_data=ld$timevar_data$death,
##'                    timevar_data=ld$timevar_data[c("bleeding","changeSBP","A","B")])
##' x <- add_baseline_data(x,data=ld$baseline_data)
##' x <- long_to_wide(x,start_followup_date=0)
##' x <- prepare_rtmle_data(x)
##' x <- protocol(x,name = "Always_A",
##'               intervention = data.frame(time=x$intervention_nodes,
##'                                         "A" = factor("1",levels = c("0","1"))))
##' x <- protocol(x,name = "Never_A",
##'               intervention = data.frame(time=x$intervention_nodes,
##'                                         "A" = factor("0",levels = c("0","1"))))
##' x <- protocol(x,name = "Use_A_not_B",
##'               intervention = data.frame(time=x$intervention_nodes,
##'                                         "A" = factor("1",levels = c("0","1")),
##'                                         "B" = factor("0",levels = c("0","1"))))
##' x <- model_formula(x)
##' x$models
##' # remove age from all formulas
##' x <- model_formula(x,exclusion_rules=list("*"="age"))
##' # remove age from the stroke_1 formula
##' x <- model_formula(x,exclusion_rules=list("stroke_1"="age"))
##' x$models
##' # remove age from all stroke formulas
##' x <- model_formula(x,exclusion_rules=list("stroke_*"="age"))
##' # remove age and changeSBP_0 from dropout_1 and dropout_2 formulas
##' x <- model_formula(x,exclusion_rules=list("dropout_[1:2]"="age|changeSBP_0"))
##' x$models
##' # remove all changeSBP_t variables from dropout_1 and dropout_2 formulas
##' x <- model_formula(x,exclusion_rules=list("dropout_[1:2]"="changeSBP_*"))
##' x$models
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
model_formula <- function(x,
                          propensity_model = "joint",
                          Markov = NULL,
                          verbose = TRUE,
                          exclude_variables = NULL,
                          exclusion_rules = NULL,
                          inclusion_rules = NULL,
                          ...){
    time_node <- NULL
    exclude_variables = c("start_followup_date",exclude_variables)
    if (length(x$protocols) == 0) {stop("No protocols registered in object, hence it is unclear which variables are intervened upon.")}
    name_time_covariates <- setdiff(x$names$name_time_covariates,exclude_variables)
    name_baseline_covariates <- setdiff(x$names$name_baseline_covariates,exclude_variables)
    if (length(name_time_covariates)>0){
        if (length(Markov)>0 && Markov[[1]]!="")
            if (any(not_found <- !(Markov%in%name_time_covariates)))
                stop(paste0("The following variables in argument Markov do not match time_covariates:\n",
                            paste(Markov[not_found],collapse=", ")))
    }
    #
    # loop across time points looking at treatments, censoring, outcomes from the beginning of the interval
    #
    model_formulas <- lapply(x$intervention_nodes,function(tk){
        # Evaluate the rule on the full prepared history. Besides making the
        # returned same-node metadata available, this preserves the callback
        # contract for rules that use cohort-level history.
        metadata_data <- data.table::as.data.table(x$prepared_data)
        protocol_specs <- lapply(names(x$protocols),function(protocol_name){
            pro <- x$protocols[[protocol_name]]
            # no intervention corresponds to setting NA or to not
            # have a line for the variable(s) in the intervention_table
            vals <- pro$intervention_table[time_node == tk][["value"]]
            if (length(vals) == 0 || all(is.na(vals))){
                NULL
            }else{
                names_vals <- pro$intervention_table[time_node == tk][["variable"]][!is.na(vals)]
                vals <- vals[!is.na(vals)]
                names(vals) <- names_vals
                intervention <- evaluate_intervention(
                    protocol = pro,
                    data = metadata_data,
                    intervention_table = pro$intervention_table,
                    time_node = tk
                )
                treatment_variables <- names(vals)

                variables <- intervention$propensity_variables
                if (length(variables) == 0) {
                    propensity_variables <- NULL
                } else if (is.character(variables)) {
                    propensity_variables <- variables
                } else {
                    joint_name <- paste(treatment_variables, collapse = ",")
                    if (joint_name %in% names(variables)) {
                        propensity_variables <- variables[[joint_name]]
                    } else if (length(treatment_variables) == 1L &&
                               ".default" %in% names(variables)) {
                        propensity_variables <- variables[[".default"]]
                    } else {
                        propensity_variables <- unique(unlist(
                            variables[intersect(treatment_variables, names(variables))],
                            use.names = FALSE
                        ))
                    }
                }

                instructions <- intervention$propensity_instructions
                if (length(instructions) == 0) {
                    propensity_instructions <- NULL
                } else {
                    joint_name <- paste(treatment_variables, collapse = ",")
                    if (joint_name %in% names(instructions)) {
                        propensity_instructions <- instructions[[joint_name]]
                    } else if (length(treatment_variables) == 1L &&
                               ".default" %in% names(instructions)) {
                        propensity_instructions <- instructions[[".default"]]
                    } else {
                        selected <- instructions[
                            intersect(treatment_variables, names(instructions))
                        ]
                        if (length(selected) == 0L) {
                            propensity_instructions <- NULL
                        } else if (length(selected) == 1L &&
                                   length(treatment_variables) == 1L) {
                            propensity_instructions <- selected[[1L]]
                        } else if (
                            length(selected) == length(treatment_variables) &&
                            all(vapply(
                                selected,
                                function(z) identical(z$mode, "adherence"),
                                logical(1)
                            ))
                        ) {
                            strata <- lapply(selected, `[[`, "stratify_by")
                            if (all(vapply(
                                strata,
                                identical,
                                logical(1),
                                strata[[1L]]
                            ))) {
                                propensity_instructions <- selected[[1L]]
                            } else {
                                propensity_instructions <- NULL
                            }
                        } else {
                            # Mixed instructions remain task-specific and are
                            # resolved during probability fitting.
                            propensity_instructions <- NULL
                        }
                    }
                }
                list(values = vals,
                     propensity_variables = propensity_variables,
                     propensity_instructions = propensity_instructions)
            }
        })
        names(protocol_specs) <- names(x$protocols)
        protocol_specs <- Filter(Negate(is.null),protocol_specs)
        all_vars <- lapply(protocol_specs,function(spec) spec$values)
        # censoring variables before outcome (no model for censoring at time zero)
        if(length(x$names$censoring)>0){
            censvalue <- paste0("'",x$names$uncensored_label,"'")
            names(censvalue) = paste0(x$names$censoring,"_",(tk+1))
            all_vars <- c(all_vars,list("censoring" = censvalue))
        }
        # outcome variables (no model for outcome at time zero)
        outvalue <- 1
        names(outvalue) <- paste0(x$names$outcome,"_",(tk+1))
        all_vars <- c(all_vars,list("outcome" = outvalue))
        # return vector of character formulas
        tk_forms <- lapply(names(all_vars), function(nav){
            vv <- all_vars[[nav]]
            ## this may be confusing but at intervention node k
            ## we evaluate treatment in the previous interval [t_{k-1},t_k] but
            ## censoring and outcome in the next interval [t_{k},t_{k+1}]
            if (nav%in% c("censoring","outcome")){
                eval_time <- tk+1
                additional_variables <- NULL
                propensity_instructions <- NULL
            } else{
                eval_time <- tk
                propensity_instructions <- protocol_specs[[nav]]$propensity_instructions
                additional_variables <- unique(c(
                    protocol_specs[[nav]]$propensity_variables,
                    if (is.null(propensity_instructions)) {
                        NULL
                    } else {
                        propensity_instructions$stratify_by
                    }
                ))
            }
            ff <- formalize(timepoint = eval_time,
                            available_names = names(x$prepared_data),
                            name_outcome_variable = names(vv),
                            outcome_value = as.character(vv),
                            name_baseline_covariates = name_baseline_covariates,
                            name_time_covariates  = name_time_covariates,
                            Markov = Markov,
                            constant_variables = x$names$name_constant_variables,
                            exclusion_rules = exclusion_rules,
                            inclusion_rules = inclusion_rules,
                            handle_concomitant_variables = propensity_model,
                            unwanted_variables = exclude_variables,
                            additional_variables = additional_variables)
            # Keep treatment task names as the model keys, while changing the
            # response for an adherence instruction to an internal response
            # column populated by intervention_probabilities().
            if (nav != "censoring" && nav != "outcome" &&
                !is.null(propensity_instructions)) {
                task_names <- names(ff)
                for (j in seq_along(ff)) {
                    ff[[j]]$propensity_mode <- propensity_instructions$mode
                    ff[[j]]$propensity_stratify_by <-
                        propensity_instructions$stratify_by
                    if (identical(propensity_instructions$mode, "adherence")) {
                        adherence_name <- paste0(
                            ".rtmle_adherence_",
                            gsub("[^A-Za-z0-9_]", "_", task_names[[j]])
                        )
                        formula_parts <- strsplit(
                            ff[[j]]$formula,
                            " ~ ",
                            fixed = TRUE
                        )[[1L]]
                        if (length(formula_parts) != 2L) {
                            stop(
                                "Could not identify the response in propensity formula: ",
                                ff[[j]]$formula
                            )
                        }
                        ff[[j]]$formula <- paste0(
                            adherence_name,
                            " ~ ",
                            formula_parts[[2L]]
                        )
                    }
                }
            }
            ff
        })
        names(tk_forms) <- names(all_vars)
        tk_forms
    })
    names(model_formulas) <- paste0("time_",x$intervention_nodes)
    x$models <- model_formulas
    x
}


######################################################################
### model_formula.R ends here
