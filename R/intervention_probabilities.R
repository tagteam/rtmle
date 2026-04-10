### intervention_probabilities.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 17 2024 (09:26) 
## Version: 
## Last-Updated: apr 10 2026 (15:40) 
##           By: Thomas Alexander Gerds
##     Update #: 519
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
                                       progressbar){
    variable = NULL
    N <- NROW(x$prepared_data)
    # restrict actions to the intervention_nodes before the max of the current run 
    action_nodes <- x$intervention_nodes[x$intervention_nodes <= max_intervention_node]
    current_protocol <- x$protocols[[protocol_name]]
    # always construct intervention_table for all intervention nodes
    intervention_table <- current_protocol$intervention_table[time <= max(x$intervention_nodes)]
    NC <- NROW(intervention_table)
    if (length(x$names$censoring)){
        NC <- NC + sum(paste0(x$names$censoring,"_",action_nodes+1) %chin% names(x$prepared_data))
    }
    #
    # define a matrix which indicates if the intervention is followed
    # the matrix should have a row for each individual and a column for
    # each intervention node (time interval)
    #
    if (
        length(intervention_match <- x$protocols[[protocol_name]]$intervention_match) == 0
    ){
        intervention_match <- matrix(0,ncol = length(x$intervention_nodes),nrow = N)
        intervention_match_names <- vector("character",length(x$intervention_nodes))
        previous <- rep(1,N)
        for(k in x$intervention_nodes){
            intervention_variables <- intervention_table[time == k][["variable"]]
            if (length(intervention_variables)>0){
                observed_values <- x$prepared_data[,intervention_variables,with = FALSE]
                for (v in 1:length(intervention_variables)){
                    # when there are multiple intervention variables
                    # all observed values must match
                    intervention_values <- intervention_table[time == k & variable == intervention_variables[[v]]][["value"]]
                    intervention_match[,k+1] <- previous <- previous*(observed_values[[intervention_variables[[v]]]] %in% intervention_values)
                    # when there are multiple treatment variables we paste-collapse the names
                    intervention_match_names[k+1] <- paste0(intervention_variables,collapse = ",")
                }
            }else{
                # FIXME: does this work when we only intervene on baseline treatment but not on subsequent treatments?
                # last value carried forward
                intervention_match[,k+1] <- previous
                intervention_match_names[k+1] <- paste0("No_intervention",k)
            }
        }
        ## the intervention_match matrix has one column per time point
        colnames(intervention_match) <- intervention_match_names
        x$protocols[[protocol_name]]$intervention_match <- intervention_match
    }
    #
    # construct the matrix with the cumulative intervention/censoring probabilities
    #
    if (refit || (NCOL(current_protocol$cumulative_intervention_probs) < NC)){
        intervention_probs <- x$prepared_data[,x$names$id,with = FALSE]
        # run backwards through time so that intervention_probs
        # has the correct column order
        intervention_last_nodes <- vector("character",length = length(action_nodes))
        if (progressbar){
            message("Fitting propensity score and censoring models: ",protocol_name)
            progress <- txtProgressBar(max = NC, style = progressbar, width=20)
            action <- 0
        }
        for (k in rev(action_nodes)){
            if (length(x$followup) == 0){
                outcome_free_and_uncensored <- rep(TRUE,N)
            }else{
                # store who is at_risk at time_k
                outcome_free_and_uncensored <- x$followup$last_interval >= k
            }
            # data used to fit the models in this time interval 
            current_data <- x$prepared_data[outcome_free_and_uncensored]
            # set the treatment variables to their protocolled values
            if (length(x$protocol[[protocol_name]]$intervene_function) == 0){
                stop(paste0("No intervene function defined for protocol ",protocol_name,"."))
            }
            intervened_data <- do.call(x$protocol[[protocol_name]]$intervene_function,
                                       list(data = current_data,
                                            intervention_table = intervention_table,time = k))
            # loop through nuisance parameter model formulas
            mk <- x$models[[paste0("time_",k)]]
            intervention_probs_k <- current_data[,x$names$id,with = FALSE]
            for (NUI in c(protocol_name,"censoring")){
                nuisance_parameter_models <- mk[[NUI]]
                for (NP in names(nuisance_parameter_models)){
                    current_formula <- nuisance_parameter_models[[NP]]$formula
                    ## print(current_formula)
                    if (progressbar){
                        action <- action + 1
                        setTxtProgressBar(progress,action)
                    }
                    nuisance_fit <- fitter(intervention_node = k,
                                           learner = learner,
                                           formula = current_formula,
                                           data = current_data,
                                           intervened_data = intervened_data,
                                           id_variable = x$names$id,
                                           minority_threshold = x$tuning_parameters$minority_threshold,
                                           seed = seed,
                                           diagnostics = x$diagnostics)
                    # store the fit
                    x$models[[paste0("time_",k)]][[NUI]][[NP]]$fit <- attr(nuisance_fit,"fit")
                    # update diagnostics
                    if (length(dia <- attr(nuisance_fit,"diagnostics"))>0){
                        if (is.null(x$diagnostics)){
                            x$diagnostics <- dia
                        }else{
                            for (dd in names(dia))
                                x$diagnostics[[dd]] <- dia[[dd]]
                        }
                    }
                    # remove attributes predicted values
                    nuisance_fit <- as.numeric(nuisance_fit)
                    # check predicted values
                    if (any(is.na(nuisance_fit))){
                        stop(paste0("Fitting nuisance parameter model returned missing values:\n",
                                    current_formula))
                    }
                    if (all(nuisance_fit == 0)){
                        stop(paste0("Nuisance parameter model is exactly zero:\n",
                                    current_formula))
                    }
                    if (any(nuisance_fit<0) || any(nuisance_fit>1)){
                        nuisance_fit <- pmax(0,pmin(1,nuisance_fit))
                        x$diagnostics$probabilities_off_range <- rbind(x$diagnostics$probabilities_off_range,
                                                                       data.table(Intervention_Node = k,
                                                                                  Nuisance_parameter = NUI,
                                                                                  Model = NP))
                    }
                    # add columns to the intervention_probs matrix
                    intervention_probs_k <- cbind(intervention_probs_k,
                                                  matrix(nuisance_fit,ncol = 1,
                                                         dimnames = list(NULL,NP)))
                }
                intervention_last_nodes[k+1] <- colnames(intervention_probs_k)[NCOL(intervention_probs_k)]
            }
            intervention_probs <- intervention_probs_k[intervention_probs,on = "id"]
        }
        # Store the intervention probabilities
        x$protocols[[protocol_name]]$intervention_probs <- intervention_probs
        x$protocols[[protocol_name]]$intervention_last_nodes <- intervention_last_nodes
        # FIXME: write this rowCumprods in armadillo
        x$protocols[[protocol_name]]$cumulative_intervention_probs <- matrixStats::rowCumprods(as.matrix(intervention_probs[,-1,with = FALSE]))
    }
    if (progressbar){cat("\n")}
    x
}

######################################################################
### intervention_probabilities.R ends here
