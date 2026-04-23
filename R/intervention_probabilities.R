### intervention_probabilities.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 17 2024 (09:26) 
## Version: 
## Last-Updated: apr 23 2026 (16:52) 
##           By: Thomas Alexander Gerds
##     Update #: 550
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
    variable = type = NULL
    # set the treatment variables to their protocolled values
    if (length(x$protocols[[protocol_name]]$intervene_function) == 0){
        stop(paste0("No intervene function defined for protocol ",protocol_name,"."))
    }
    N <- NROW(x$prepared_data)
    # restrict actions to the intervention_nodes before the max of the current run
    action_nodes <- x$intervention_nodes[x$intervention_nodes <= max_intervention_node]
    current_protocol <- x$protocols[[protocol_name]]
    # extract intervention_table for the intervention nodes before max_intervention_node
    intervention_table <- current_protocol$intervention_table[time <= max_intervention_node]
    # the number of columns is defined by the number of rows in the intervention table
    NC <- NROW(intervention_table)
    # plus the number of censoring variables
    if (length(x$names$censoring)){
        NC <- NC + sum(paste0(x$names$censoring,"_",action_nodes+1) %chin% names(x$prepared_data))
    }
    #
    # construct a matrices with the intervention/censoring probabilities
    #
    if (refit ||
        # only run the necessary models for the current maximal time
        # horizon which is here defined by NC via max_intervention_node 
        (NCOL(current_protocol$cumulative_intervention_probs) < NC)){
        if (progressbar){
            message("Fitting propensity score and censoring models: ",protocol_name)
            progress <- txtProgressBar(max = NC, style = progressbar, width=20)
            action <- 0
        }
        intervention_probs <- matrix(NA_real_,nrow = N,ncol = NC)
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
        colnames(intervention_probs) <- task_list$variable
        intervention_last_nodes <- task_list[,variable[.N],by = time]$V1
        task_column <- NC
        for (task in 1:nrow(task_list)){
            k <- task_list[task,time]
            if (length(x$followup) == 0){
                outcome_free_and_uncensored <- rep(TRUE,N)
            }else{
                # store who is at_risk at time_k
                outcome_free_and_uncensored <- x$followup$last_interval >= k
            }
            # prepare data used to fit the models in this time interval 
            current_data <- x$prepared_data[outcome_free_and_uncensored]
            # prepare data used to predict the intervention/censoring probabilities
            intervened_data <- do.call(x$protocols[[protocol_name]]$intervene_function,
                                       list(data = current_data,
                                            intervention_table = intervention_table,
                                            time = k))
            # fit all nuisance parameter models for intervention node k (time interval k)
            ## current_formula <- x$models[[paste0("time_",k)]][[task_list[task,type]]][[task_list[task,variable]]]$formula
            current_formula <- task_list[task,formula]
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
            x$models[[paste0("time_",k)]][[task_list[task,type]]][[task_list[task,variable]]]$fit <- attr(nuisance_fit,"fit")
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
                                                               task_list[k])
            }
            # add columns to the intervention_probs matrix
            intervention_probs[outcome_free_and_uncensored,task] <- nuisance_fit
        }
        # Store the intervention probabilities
        x$protocols[[protocol_name]]$intervention_probs <- intervention_probs
        x$protocols[[protocol_name]]$intervention_last_nodes <- intervention_last_nodes
        # FIXME: write this rowCumprods in armadillo
        #        and only keep the columns of the intervention_last_nodes
        x$protocols[[protocol_name]]$cumulative_intervention_probs <- matrixStats::rowCumprods(as.matrix(intervention_probs))
    }
    if (progressbar){cat("\n")}
    x
}

######################################################################
### intervention_probabilities.R ends here
