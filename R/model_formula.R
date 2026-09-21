### model_formula.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 16 2025 (08:58) 
## Version: 
## Last-Updated: sep 19 2026 (07:34)
##           By: Thomas Alexander Gerds
##     Update #: 168
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
##' \code{L_k}. Dynamic adherence models therefore use the history available
##' before the treatment decision. See \code{adherence_model_strata} in
##' \code{\link{regime}} for optional subgroup-specific fitting; those
##' variables are not added as current-node predictors.
##'
##' The treatment factorization is set separately for each regime with
##' \code{multiple_treatment_factorization}. See \code{\link{regime}}
##' for the dynamic-regime instruction interface.
##' The factorization choice matters only when multiple treatment variables
##' are assigned at the same intervention node.
##'
##' For a dynamic regime, a user-supplied \code{intervene_function} changes the
##' treatment response to observed-versus-intervened adherence. The response is
##' binary, and \code{adherence_model_strata} in \code{\link{regime}} can
##' request separate fits in pre-decision subgroups. Stratification variables
##' are not added as current-interval predictors.
##'
##' @title Model formulas for nuisance parameters
##' @param x An object of class \code{"rtmle"} with prepared data and at least
##'   one regime defined by \code{\link{regime}}.
##' @param Markov Names of time-dependent variables that should appear only with
##'   their most recent values on the right-hand side of formulas.
##' @param exclude_variables Variables to exclude from all formulas of the nuisance parameter models.
##' @param exclusion_rules Additional exclusion rules given as a
##'   named list. Names are variables occurring on the left-hand side of a
##'   formula, and elements are variables to exclude from the right-hand side.
##'   E.g., \code{list("A_4"="age")} removes variable \code{"age"} from the propensity
##'   score model for treatment variable \code{"A_4"}.
##'   Regular expressions are interpreted by \code{grep} so that
##'   \code{list("*"="age")} removes variable \code{"age"} from all formulas and
##'   \code{list("dropout_*"="age")} removes variable \code{"age"} from all formulas
##'   where the outcome is a censoring node (assuming \code{x$name$censoring="dropout"}).
##' @param inclusion_rules Additional inclusion rules given as a
##'   named list. Names are variables occurring on the left-hand side of a
##'   formula, and elements are variables to include on the right-hand side.
##' @param verbose Logical. If \code{FALSE} suppress all messages. \code{TRUE} is the default.
##' @return The modified \code{rtmle} object.
##' @seealso \code{\link{prepare_rtmle_data}}, \code{\link{regime}},
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
##' x <- discretize_data(x,start_followup_date=0)
##' x <- prepare_rtmle_data(x)
##' x <- regime(x,name = "Always_A",
##'               intervention = data.frame(time=x$intervention_nodes,
##'                                         "A" = factor("1",levels = c("0","1"))))
##' x <- regime(x,name = "Never_A",
##'               intervention = data.frame(time=x$intervention_nodes,
##'                                         "A" = factor("0",levels = c("0","1"))))
##' x <- regime(x,name = "Use_A_not_B",
##'               intervention = data.frame(time=x$intervention_nodes,
##'                                         "A" = factor("1",levels = c("0","1")),
##'                                         "B" = factor("0",levels = c("0","1"))),
##'               multiple_treatment_factorization = "sequential")
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
                          Markov = NULL,
                          verbose = TRUE,
                          exclude_variables = NULL,
                          exclusion_rules = NULL,
                          inclusion_rules = NULL){
    time_node <- NULL
    default_propensity_factorization <- "joint"
    exclude_variables = c("start_followup_date",exclude_variables)
    if (length(x$regimes) == 0) {stop("No regimes registered in object, hence it is unclear which variables are intervened upon.")}
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
        # Evaluate the rule on the full prepared history so that the
        # intervention-updated treatment values are resolved consistently
        # while formulas are built.
        metadata_data <- data.table::as.data.table(x$prepared_data)
        regime_specs <- lapply(names(x$regimes),function(regime_name){
            pro <- x$regimes[[regime_name]]
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
                    regime = pro,
                    data = metadata_data,
                    intervention_table = pro$intervention_table,
                    time_node = tk
                )
                treatment_variables <- names(vals)

                adherence_model_strata <-
                    intervention$adherence_model_strata
                dynamic_adherence <- intervention$dynamic_adherence
                regime_factorization <-
                    intervention$multiple_treatment_factorization
                validate_multiple_treatment_factorization(
                    regime_factorization,
                    paste0(
                        "`multiple_treatment_factorization` for regime '",
                        regime_name, "'"
                    )
                )
                list(values = vals,
                     adherence_model_strata = adherence_model_strata,
                     dynamic_adherence = dynamic_adherence,
                     propensity_factorization = regime_factorization)
            }
        })
        names(regime_specs) <- names(x$regimes)
        regime_specs <- Filter(Negate(is.null),regime_specs)
        all_vars <- lapply(regime_specs,function(spec) spec$values)
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
                adherence_model_strata <- NULL
                dynamic_adherence <- FALSE
                task_propensity_factorization <-
                    default_propensity_factorization
            } else{
                eval_time <- tk
                task_propensity_factorization <-
                    regime_specs[[nav]]$propensity_factorization
                adherence_model_strata <-
                    regime_specs[[nav]]$adherence_model_strata
                dynamic_adherence <- regime_specs[[nav]]$dynamic_adherence
                additional_variables <- NULL
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
                            handle_concomitant_variables =
                                task_propensity_factorization,
                            unwanted_variables = exclude_variables,
                            additional_variables = additional_variables)
            # Keep treatment task names as the model keys, while changing the
            # response for a dynamic intervention to an internal adherence
            # column populated by intervention_probabilities().
            if (nav != "censoring" && nav != "outcome" &&
                isTRUE(dynamic_adherence)) {
                model_is_flat <- length(ff) == 0L ||
                    all(vapply(
                        ff,
                        function(model) !is.null(model$formula),
                        logical(1)
                    ))
                if (model_is_flat) {
                    model_names <- names(ff)
                    for (j in seq_along(ff)) {
                        ff[[j]]$dynamic_adherence <- TRUE
                        ff[[j]]$adherence_model_strata <-
                            adherence_model_strata
                        {
                            adherence_name <- paste0(
                                ".rtmle_adherence_",
                                gsub(
                                    "[^A-Za-z0-9_]",
                                    "_",
                                    model_names[[j]]
                                )
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
                } else {
                    # Sequential propensity models are nested by ordered
                    # treatment step, so the task name is held by each inner
                    # list rather than by ff itself.
                    for (j in seq_along(ff)) {
                        model_names <- names(ff[[j]])
                        for (h in seq_along(ff[[j]])) {
                            ff[[j]][[h]]$dynamic_adherence <- TRUE
                            ff[[j]][[h]]$adherence_model_strata <-
                                adherence_model_strata
                            {
                                adherence_name <- paste0(
                                    ".rtmle_adherence_",
                                    gsub(
                                        "[^A-Za-z0-9_]",
                                        "_",
                                        model_names[[h]]
                                    )
                                )
                                formula_parts <- strsplit(
                                    ff[[j]][[h]]$formula,
                                    " ~ ",
                                    fixed = TRUE
                                )[[1L]]
                                if (length(formula_parts) != 2L) {
                                    stop(
                                        "Could not identify the response in propensity formula: ",
                                        ff[[j]][[h]]$formula
                                    )
                                }
                                ff[[j]][[h]]$formula <- paste0(
                                    adherence_name,
                                    " ~ ",
                                    formula_parts[[2L]]
                                )
                            }
                        }
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
