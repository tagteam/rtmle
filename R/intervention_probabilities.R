### intervention_probabilities.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 17 2024 (09:26) 
## Version: 
## Last-Updated: jun 18 2026 (08:52) 
##           By: Thomas Alexander Gerds
##     Update #: 613
#----------------------------------------------------------------------
## 
### Commentary: 
#
# fit all treatment and censoring nuisance parameter models
# and gather propensitity scores and censoring probabilities in a matrix
# such that censoring comes last at each intervention node
#
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
intervention_probabilities <- function(x,
                                       regime_name,
                                       max_intervention_node,
                                       refit = FALSE,
                                       learner,
                                       seed,
                                       progressbar,
                                       save_fitted_objects = FALSE){
    variable = type = time_node = value = NULL
    # set the treatment variables to their regimeled values
    if (length(x$regimes[[regime_name]]$intervene_function) == 0){
        stop(paste0("No intervene function defined for regime ",regime_name,"."))
    }
    N <- NROW(x$prepared_data)
    # restrict actions to the intervention_nodes before the max of the current run
    action_nodes <- x$intervention_nodes[x$intervention_nodes <= max_intervention_node]
    current_regime <- x$regimes[[regime_name]]
    # extract intervention_table for the intervention nodes before max_intervention_node
    intervention_table <- na.omit(current_regime$intervention_table[time_node <= max_intervention_node])
    #
    # construct a matrices with the intervention/censoring probabilities
    #
    task_list <- do.call(rbind,lapply(action_nodes, function(k){
        regime_models <- x$models[[paste0("time_", k)]][[regime_name]]
        treatment_tasks <- if (length(regime_models) == 0L) {
            NULL
        } else if (all(vapply(
            regime_models,
            function(model) !is.null(model$formula),
            logical(1)
        ))) {
            # Joint and independent propensity models are represented as a
            # named list of model specifications.
            do.call(rbind, lapply(names(regime_models), function(v) {
                data.table(
                    time = k,
                    type = regime_name,
                    variable = v,
                    formula = regime_models[[v]]$formula
                )
            }))
        } else {
            # Sequential propensity models are represented as a list of
            # ordered steps, each containing one named model specification.
            do.call(rbind, lapply(regime_models, function(step) {
                do.call(rbind, lapply(names(step), function(v) {
                    data.table(
                        time = k,
                        type = regime_name,
                        variable = v,
                        formula = step[[v]]$formula
                    )
                }))
            }))
        }
        censoring_tasks <- do.call(rbind, lapply(
            x$models[[paste0("time_", k)]]["censoring"],
            function(w) {
                do.call(rbind, lapply(names(w), function(v) {
                    data.table(
                        time = k,
                        type = "censoring",
                        variable = v,
                        formula = w[[v]]$formula
                    )
                }))
            }
        ))
        rbind(treatment_tasks, censoring_tasks)
    }))
    # the number of columns is defined by the number of censoring models plus the
    # number of propensitity scores models which in case of multiple treatment variables
    # depends on the type of propensity score modelling (joint vs sequential)
    NC <- NROW(task_list)
    if (refit ||
        # only run the necessary models for the current maximal time
        # horizon which is here defined by NC via max_intervention_node 
        (NCOL(current_regime$cumulative_intervention_probs) < NC)){
        if (progressbar){
            message("Fitting propensity score and censoring models: ",regime_name)
            progress <- txtProgressBar(max = NC, style = progressbar, width=20)
            action <- 0
        }
        # FIXME: if the first time points were readily run but not yet all then
        #        (unless refit is TRUE) we would like to preserve the probs
        intervention_probs <- matrix(NA_real_,nrow = N,ncol = NC)
        colnames(intervention_probs) <- task_list$variable
        # this following code is slightly cryptic and should be improved.
        # the aim is to make the functionality robust in situations where there is no weighting
        # or no intervention at some time nodes
        ipw_last_nodes <- rep(NA,length(action_nodes))
        names(ipw_last_nodes) <- paste0("node_",action_nodes)
        ipw_last_nodes_data <- task_list[,variable[.N],by = time]
        if (NROW(ipw_last_nodes_data)>0){
            actual_ipw_last_nodes <- ipw_last_nodes_data$V1
            names(actual_ipw_last_nodes) <- paste0("node_",ipw_last_nodes_data$time)
            ipw_last_nodes[names(actual_ipw_last_nodes)] <- actual_ipw_last_nodes
        }
        # now the same for non-censoring intervention nodes
        intervention_last_nodes <- rep(NA,length(action_nodes))
        names(intervention_last_nodes) <- paste0("node_",action_nodes)
        # intervention_match has one column per intervention node containing
        # all treatment variables separated by commas (for example
        # A_1,B_1).  This differs from the last nuisance task under
        # sequential or independent propensity models, where the last task
        # would be only B_1.  Keep the matching name aligned with the
        # intervention table so the TMLE update can find the column.
        intervention_last_nodes_data <- intervention_table[
            !is.na(value),
            list(V1 = paste(variable, collapse = ",")),
            by = time_node
        ]
        data.table::setnames(
            intervention_last_nodes_data,
            "time_node",
            "time"
        )
        if (NROW(intervention_last_nodes_data)>0){
            actual_intervention_last_nodes <- intervention_last_nodes_data$V1
            names(actual_intervention_last_nodes) <- paste0("node_",intervention_last_nodes_data$time)
            intervention_last_nodes[names(actual_intervention_last_nodes)] <- actual_intervention_last_nodes
        }
        evaluated_interventions <- list()
        for (task in seq_len(nrow(task_list))){
            k <- task_list[task,time]
            if (length(x$followup) == 0){
                outcome_free_and_uncensored <- rep(TRUE,N)
            }else{
                # store who is at_risk at time_k
                outcome_free_and_uncensored <- x$followup$last_interval >= k
            }
            # prepare data used to fit the models in this time interval 
            current_data <- x$prepared_data[outcome_free_and_uncensored]
            # Evaluate a user-defined intervention once per node. A custom
            # intervention function requests an adherence model; optional
            # adherence_model_strata only controls whether that model is fit
            # separately in pre-decision subgroups.
            intervention_key <- paste0("node_",k)
            if (is.null(evaluated_interventions[[intervention_key]])){
                evaluated_interventions[[intervention_key]] <- evaluate_intervention(
                    regime = current_regime,
                    data = current_data,
                    intervention_table = intervention_table,
                    time_node = k,
                    full_n = N,
                    row_indices = which(outcome_free_and_uncensored)
                )
            }
            intervention <- evaluated_interventions[[intervention_key]]
            intervened_data <- data.table::copy(intervention$data)
            task_variable <- as.character(task_list[task,variable])
            dynamic_adherence <- FALSE
            adherence_model_strata <- NULL
            current_formula <- as.character(task_list[task,formula])
            if (task_list[task,type] != "censoring"){
                current_treatment_variables <- intervention_table[
                    time_node == k
                ][["variable"]]
                task_treatment_variables <- intersect(
                    strsplit(task_variable, ",", fixed = TRUE)[[1L]],
                    current_treatment_variables
                )
                dynamic_adherence <- isTRUE(intervention$dynamic_adherence)
                adherence_model_strata <- intervention$adherence_model_strata
                if (dynamic_adherence) {
                    task_variables <- strsplit(
                        task_variable,
                        ",",
                        fixed = TRUE
                    )[[1L]]
                    adherence <- rep(TRUE, NROW(current_data))
                    for (variable_name in task_variables) {
                        if (!(variable_name %in% names(current_data)) ||
                            !(variable_name %in% names(intervened_data))) {
                            stop(
                                "Cannot construct adherence for treatment variable ",
                                variable_name,
                                "."
                            )
                        }
                        matches <- as.character(current_data[[variable_name]]) ==
                            as.character(intervened_data[[variable_name]])
                        matches[is.na(matches)] <- FALSE
                        adherence <- adherence & matches
                    }
                    adherence <- as.integer(adherence)
                    adherence_name <- paste0(
                        ".rtmle_adherence_",
                        gsub("[^A-Za-z0-9_]", "_", task_variable)
                    )
                    current_data[[adherence_name]] <- adherence
                    # The response is not needed for prediction, but retaining
                    # it makes the newdata contract safe for all learners.
                    intervened_data[[adherence_name]] <- 1L
                    formula_parts <- strsplit(
                        current_formula,
                        " ~ ",
                        fixed = TRUE
                    )[[1L]]
                    if (length(formula_parts) != 2L) {
                        stop("Could not identify the response in propensity formula: ",
                             current_formula)
                    }
                    current_formula <- paste0(
                        adherence_name,
                        " ~ ",
                        formula_parts[[2L]]
                    )
                } else {
                    # The ordinary treatment formula has a nominal, fixed
                    # response. The default intervention function is static;
                    # a custom dynamic intervention is handled above.
                    changed_from_nominal <- rep(FALSE, NROW(current_data))
                    for (treatment_variable in task_treatment_variables){
                        nominal_value <- intervention_table[
                            time_node == k & variable == treatment_variable
                        ][["value"]]
                        intervened_value <- intervention$data[[treatment_variable]]
                        changed <- as.character(intervened_value) !=
                            as.character(nominal_value[[1L]])
                        changed[is.na(changed)] <- TRUE
                        changed_from_nominal <- changed_from_nominal | changed
                    }
                    if (any(changed_from_nominal)){
                        stop(
                            "The intervention changes a nominal treatment value for ",
                            task_variable,
                            " without a custom `intervene_function`."
                        )
                    }
                }
            }
            # Fit all nuisance parameter models for intervention node k. A
            # dynamic adherence model can be split by pre-decision strata.
            if (task_list[task,type] == "censoring" ||
                !dynamic_adherence ||
                is.null(adherence_model_strata) ||
                length(adherence_model_strata) == 0L) {
                fit_groups <- list(.all = seq_len(NROW(current_data)))
            } else {
                strata <- adherence_model_strata
                if (length(strata) == 0L) {
                    fit_groups <- list(.all = seq_len(NROW(current_data)))
                } else {
                    missing_strata <- setdiff(strata, names(current_data))
                    if (length(missing_strata) > 0L) {
                        stop(
                            "Unknown `strata` variable(s): ",
                            paste(missing_strata, collapse = ", "),
                            "."
                        )
                    }
                    if (any(vapply(
                        strata,
                        function(variable_name) anyNA(current_data[[variable_name]]),
                        logical(1)
                    ))) {
                        stop(
                            "`strata` variables must be observed for every row used ",
                            "to fit an adherence propensity."
                        )
                    }
                    strata_values <- lapply(
                        strata,
                        function(variable_name) {
                            as.character(current_data[[variable_name]])
                        }
                    )
                    strata_key <- do.call(
                        paste,
                        c(strata_values, sep = "\r")
                    )
                    fit_groups <- split(
                        seq_len(NROW(current_data)),
                        strata_key,
                        drop = TRUE
                    )
                }
                fit_groups <- lapply(fit_groups, identity)
            }
            if (progressbar){
                action <- action + 1
                setTxtProgressBar(progress,action)
            }
            # save censoring model
            reuse_fit <- NULL
            save_current_fit <- save_fitted_objects
            if (task_list[task,type] == "censoring"){
                if (regime_name == x$censoring_use_regime){
                    # store the fit
                    save_current_fit <- TRUE
                }else{
                    # reuse the fit
                    save_current_fit <- FALSE
                    reuse_fit <- x$models[[paste0("time_",k)]][[task_list[task,type]]][[task_list[task,variable]]][c("fit","fit_summary")]
                }
            }
            predicted_values <- rep(NA_real_, NROW(current_data))
            fitted_objects <- list()
            fit_summaries <- list()
            fit_diagnostics <- list()
            for (group_name in names(fit_groups)) {
                group_rows <- fit_groups[[group_name]]
                if (length(group_rows) == 0L) next
                nuisance_fit <- fitter(
                    intervention_node = k,
                    learner = learner,
                    formula = current_formula,
                    data = current_data[group_rows],
                    intervened_data = intervened_data[group_rows],
                    id_variable = x$names$id,
                    minority_threshold = x$tuning_parameters$minority_threshold,
                    seed = seed,
                    diagnostics = x$diagnostics,
                    save_fitted_objects = save_current_fit,
                    reuse_fit = reuse_fit
                )
                predicted_values[group_rows] <- nuisance_fit$predicted_values
                # Keep a named element even when a learner does not return a
                # fit summary (for example learn_xgboost()). Assigning with
                # [[ ]] would drop a NULL and make a successful fit look as if
                # no model had been fitted.
                fit_summaries[group_name] <- list(nuisance_fit$fit_summary)
                if (save_current_fit) {
                    fitted_objects[[group_name]] <- nuisance_fit$fit
                }
                if (length(nuisance_fit$diagnostics) > 0) {
                    fit_diagnostics[[group_name]] <- nuisance_fit$diagnostics
                }
            }
            fitted <- length(fit_summaries) > 0L
            if (!fitted) {
                fit_summary <- structure(
                    "No model fitted for this treatment task.",
                    class = "no_propensity_model"
                )
                fitted_objects <- NULL
            } else if (length(fit_summaries) == 1L &&
                       ".all" %in% names(fit_summaries)) {
                fit_summary <- fit_summaries[[".all"]]
                if (!save_current_fit) fitted_objects <- NULL
            } else {
                fit_summary <- structure(
                    fit_summaries,
                    class = c("stratified_adherence_probability", "list")
                )
                if (!save_current_fit) fitted_objects <- NULL
            }
            # store the fit
            model_time <- paste0("time_", k)
            model_type <- as.character(task_list[task, type])
            model_container <- x$models[[model_time]][[model_type]]
            model_is_flat <- length(model_container) == 0L ||
                all(vapply(
                    model_container,
                    function(model) !is.null(model$formula),
                    logical(1)
                ))
            if (model_is_flat) {
                model_index <- task_variable
            } else {
                model_index <- which(vapply(
                    model_container,
                    function(step) task_variable %in% names(step),
                    logical(1)
                ))
                if (length(model_index) != 1L) {
                    stop(
                        "Could not locate the sequential propensity model for ",
                        task_variable, " at node ", k, "."
                    )
                }
            }
            if (save_current_fit){
                fit_to_store <- fitted_objects
                if (length(fitted_objects) == 1L &&
                    ".all" %in% names(fitted_objects)) {
                    fit_to_store <- fitted_objects[[".all"]]
                }
                if (model_is_flat) {
                    model_container[[model_index]]$fit <- fit_to_store
                } else {
                    model_container[[model_index]][[task_variable]]$fit <-
                        fit_to_store
                }
            }
            if (model_is_flat) {
                model_container[[model_index]]$fit_summary <- fit_summary
            } else {
                model_container[[model_index]][[task_variable]]$fit_summary <-
                    fit_summary
            }
            x$models[[model_time]][[model_type]] <- model_container
            # update diagnostics
            for (dia in fit_diagnostics) {
                if (is.null(x$diagnostics)){
                    x$diagnostics <- dia
                }else{
                    for (dd in names(dia)){
                        x$diagnostics[[dd]] <- dia[[dd]]
                    }
                }
            }
            # check predicted values
            if (any(is.na(predicted_values))){
                stop(paste0("Fitting nuisance parameter model returned missing values:\n",
                            current_formula))
            }
            if (fitted && length(predicted_values) > 0 &&
                all(predicted_values == 0)){
                stop(paste0("Nuisance parameter model predictions are exactly zero:\n",
                            current_formula))
            }
            if (any(predicted_values<0) || any(predicted_values>1)){
                predicted_values <- pmax(0,pmin(1,predicted_values))
                x$diagnostics$probabilities_off_range <- rbind(x$diagnostics$probabilities_off_range,
                                                               task_list[task])
            }
            # add columns to the intervention_probs matrix
            intervention_probs[outcome_free_and_uncensored,task] <- predicted_values
        }
        # Store the intervention probabilities
        x$regimes[[regime_name]]$intervention_probs <- intervention_probs
        x$regimes[[regime_name]]$ipw_last_nodes <- ipw_last_nodes
        x$regimes[[regime_name]]$intervention_last_nodes <- intervention_last_nodes
        # FIXME: write this rowCumprods in armadillo
        #        and only keep the columns of the ipw_last_nodes
        x$regimes[[regime_name]]$cumulative_intervention_probs <- matrixStats::rowCumprods(as.matrix(intervention_probs))
    }
    if (progressbar){cat("\n")}
    x
}

######################################################################
### intervention_probabilities.R ends here
