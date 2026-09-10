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
                                       protocol_name,
                                       max_intervention_node,
                                       refit = FALSE,
                                       learner,
                                       seed,
                                       progressbar,
                                       save_fitted_objects = FALSE){
    variable = type = time_node = NULL
    # set the treatment variables to their protocolled values
    if (length(x$protocols[[protocol_name]]$intervene_function) == 0){
        stop(paste0("No intervene function defined for protocol ",protocol_name,"."))
    }
    N <- NROW(x$prepared_data)
    # restrict actions to the intervention_nodes before the max of the current run
    action_nodes <- x$intervention_nodes[x$intervention_nodes <= max_intervention_node]
    current_protocol <- x$protocols[[protocol_name]]
    # extract intervention_table for the intervention nodes before max_intervention_node
    intervention_table <- na.omit(current_protocol$intervention_table[time_node <= max_intervention_node])
    #
    # construct a matrices with the intervention/censoring probabilities
    #
    task_list <- do.call(rbind,lapply(action_nodes, function(k){
        rbind(do.call(rbind,lapply(x$models[[paste0("time_",k)]][protocol_name],function(w){
            do.call(rbind,lapply(names(w),function(v){
                data.table(time = k,type = protocol_name,variable = v,formula = w[[v]]$formula)
            }))})),
            do.call(rbind,lapply(x$models[[paste0("time_",k)]]["censoring"],function(w){
                do.call(rbind,lapply(names(w),function(v){
                    data.table(time = k,type = "censoring",variable = v,formula = w[[v]]$formula)
                }))})))
    }))
    # the number of columns is defined by the number of censoring models plus the
    # number of propensitity scores models which in case of multiple treatment variables
    # depends on the type of propensity score modelling (joint vs sequential)
    NC <- NROW(task_list)
    if (refit ||
        # only run the necessary models for the current maximal time
        # horizon which is here defined by NC via max_intervention_node 
        (NCOL(current_protocol$cumulative_intervention_probs) < NC)){
        if (progressbar){
            message("Fitting propensity score and censoring models: ",protocol_name)
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
        intervention_last_nodes_data <- task_list[type != "censoring",variable[.N],by = time]
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
            # Evaluate a user-defined intervention once per node. The returned
            # object may provide fixed probabilities or an adherence-model
            # instruction for the treatment task.
            intervention_key <- paste0("node_",k)
            if (is.null(evaluated_interventions[[intervention_key]])){
                evaluated_interventions[[intervention_key]] <- evaluate_intervention(
                    protocol = current_protocol,
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
            propensity_instruction <- NULL
            fixed_probability <- rep(NA_real_, NROW(current_data))
            current_formula <- as.character(task_list[task,formula])
            if (task_list[task,type] != "censoring"){
                current_treatment_variables <- intervention_table[
                    time_node == k
                ][["variable"]]
                task_treatment_variables <- intersect(
                    strsplit(task_variable, ",", fixed = TRUE)[[1L]],
                    current_treatment_variables
                )
                instructions <- intervention$propensity_instructions
                if (length(instructions) > 0L) {
                    task_variables <- strsplit(
                        task_variable,
                        ",",
                        fixed = TRUE
                    )[[1L]]
                    if (task_variable %in% names(instructions)) {
                        propensity_instruction <- instructions[[task_variable]]
                    } else if (length(task_variables) == 1L &&
                               ".default" %in% names(instructions)) {
                        propensity_instruction <- instructions[[".default"]]
                    } else {
                        selected <- instructions[
                            intersect(task_variables, names(instructions))
                        ]
                        if (length(selected) == 1L &&
                            length(task_variables) == 1L) {
                            propensity_instruction <- selected[[1L]]
                        } else if (length(selected) == 0L) {
                            propensity_instruction <- NULL
                        } else {
                            if (length(selected) != length(task_variables)) {
                                stop(
                                    "A joint propensity instruction must name the joint task or every treatment in it."
                                )
                            }
                            modes <- vapply(
                                selected,
                                `[[`,
                                character(1),
                                "mode"
                            )
                            if (all(modes == "adherence")) {
                                strata <- lapply(selected, `[[`, "stratify_by")
                                if (!all(vapply(
                                    strata,
                                    identical,
                                    logical(1),
                                    strata[[1L]]
                                ))) {
                                    stop(
                                        "Joint adherence instructions must use the same `stratify_by` variables for every treatment."
                                    )
                                }
                                propensity_instruction <- list(
                                    mode = "adherence",
                                    probability = NULL,
                                    stratify_by = strata[[1L]]
                                )
                            } else if (all(modes == "fixed")) {
                                values <- as.data.frame(
                                    lapply(selected, `[[`, "probability"),
                                    check.names = FALSE
                                )
                                overridden <- !is.na(as.matrix(values))
                                partial <- rowSums(overridden) > 0 &
                                    rowSums(overridden) < NCOL(values)
                                if (any(partial)) {
                                    stop(
                                        "Partial row-wise fixed probabilities are not valid for a joint propensity model; use sequential or independent propensity models."
                                    )
                                }
                                combined <- rep(NA_real_, NROW(values))
                                complete <- rowSums(overridden) == NCOL(values)
                                if (any(complete)) {
                                    complete_values <- as.matrix(
                                        values[complete, , drop = FALSE]
                                    )
                                    if (any(complete_values != 1)) {
                                        stop(
                                            "Non-unit joint fixed probabilities must be supplied in a column named for the joint propensity task."
                                        )
                                    }
                                    combined[complete] <- 1
                                }
                                propensity_instruction <- list(
                                    mode = "fixed",
                                    probability = combined,
                                    stratify_by = NULL
                                )
                            } else {
                                stop(
                                    "Joint propensity instructions cannot mix `fixed` and `adherence` modes."
                                )
                            }
                        }
                    }
                }
                if (!is.null(propensity_instruction) &&
                    identical(propensity_instruction$mode, "fixed")) {
                    fixed_probability <- propensity_instruction$probability
                    if (length(fixed_probability) == 1L &&
                        NROW(current_data) != 1L) {
                        fixed_probability <- rep(
                            fixed_probability,
                            NROW(current_data)
                        )
                    }
                    if (length(fixed_probability) != NROW(current_data)) {
                        stop(
                            "A fixed propensity instruction must have length 1 or the number of rows being evaluated."
                        )
                    }
                    fixed_probability <- as.numeric(fixed_probability)
                }
                if (!is.null(propensity_instruction) &&
                    identical(propensity_instruction$mode, "adherence")) {
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
                    # response (for example I(A_k == 1)). A dynamic rule that
                    # changes that value must therefore either provide a fixed
                    # instruction or explicitly request an adherence model.
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
                    if (any(changed_from_nominal & is.na(fixed_probability))){
                        stop(
                            "The intervention changes a nominal treatment value for ",
                            task_variable,
                            " without a propensity instruction. Return a fixed ",
                            "probability for those rows or use `mode = \"adherence\"`."
                        )
                    }
                }
            }
            fit_rows <- is.na(fixed_probability)
            # Fit all nuisance parameter models for intervention node k (time
            # interval k). Adherence instructions can optionally split the
            # fitting rows into separate strata.
            if (task_list[task,type] == "censoring" ||
                is.null(propensity_instruction) ||
                !identical(propensity_instruction$mode, "adherence")) {
                fit_groups <- list(.all = which(fit_rows))
            } else {
                stratify_by <- propensity_instruction$stratify_by
                if (length(stratify_by) == 0L) {
                    fit_groups <- list(.all = seq_len(NROW(current_data)))
                } else {
                    missing_strata <- setdiff(stratify_by, names(current_data))
                    if (length(missing_strata) > 0L) {
                        stop(
                            "Unknown `stratify_by` variable(s): ",
                            paste(missing_strata, collapse = ", "),
                            "."
                        )
                    }
                    if (any(vapply(
                        stratify_by,
                        function(variable_name) anyNA(current_data[[variable_name]]),
                        logical(1)
                    ))) {
                        stop(
                            "`stratify_by` variables must be observed for every row used ",
                            "to fit an adherence propensity."
                        )
                    }
                    strata_values <- lapply(
                        stratify_by,
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
                fit_groups <- lapply(fit_groups, function(indices) {
                    intersect(indices, which(fit_rows))
                })
            }
            if (progressbar){
                action <- action + 1
                setTxtProgressBar(progress,action)
            }
            # save censoring model
            reuse_fit <- NULL
            save_current_fit <- save_fitted_objects
            if (task_list[task,type] == "censoring"){
                if (protocol_name == x$censoring_use_protocol){
                    # store the fit
                    save_current_fit <- TRUE
                }else{
                    # reuse the fit
                    save_current_fit <- FALSE
                    reuse_fit <- x$models[[paste0("time_",k)]][[task_list[task,type]]][[task_list[task,variable]]][c("fit","fit_summary")]
                }
            }
            predicted_values <- fixed_probability
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
                    "No model fitted: all probabilities supplied by the propensity instruction.",
                    class = "fixed_probability_instruction"
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
            if (save_current_fit){
                fit_to_store <- fitted_objects
                if (length(fitted_objects) == 1L &&
                    ".all" %in% names(fitted_objects)) {
                    fit_to_store <- fitted_objects[[".all"]]
                }
                x$models[[paste0("time_",k)]][[task_list[task,type]]][[task_variable]]$fit <- fit_to_store
            }
            x$models[[paste0("time_",k)]][[task_list[task,type]]][[task_variable]]$fit_summary <- fit_summary
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
        x$protocols[[protocol_name]]$intervention_probs <- intervention_probs
        x$protocols[[protocol_name]]$ipw_last_nodes <- ipw_last_nodes
        x$protocols[[protocol_name]]$intervention_last_nodes <- intervention_last_nodes
        # FIXME: write this rowCumprods in armadillo
        #        and only keep the columns of the ipw_last_nodes
        x$protocols[[protocol_name]]$cumulative_intervention_probs <- matrixStats::rowCumprods(as.matrix(intervention_probs))
    }
    if (progressbar){cat("\n")}
    x
}

######################################################################
### intervention_probabilities.R ends here
