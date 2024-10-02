### run_rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  1 2024 (09:11) 
## Version: 
## Last-Updated: Oct  2 2024 (09:48) 
##           By: Thomas Alexander Gerds
##     Update #: 308
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @export
run_rtmle <- function(x,
                      targets,
                      learner = "learn_glm",
                      time_horizon,
                      refit = FALSE,
                      ...){
    requireNamespace("riskRegression")
    # for loop across targets
    available_targets <- names(x$targets)
    if (!missing(targets)) {
        if (!(all(targets %in% available_targets)))
            stop(paste0("Requested targets: \n",paste(targets,collapse = ","),"\n\navailable targets:\n",paste(available_targets,collapse = ",")))
        run_these_targets <- intersect(targets,available_targets)
    }else{
        run_these_targets <- available_targets
    }
    if (missing(time_horizon)) {
        time_horizon <- max(x$time)
    } else {
        stopifnot(all(time_horizon <= max(x$time) & time_horizon>0))
    }
    name_time_horizon <- paste0("time_horizon_",time_horizon)
    x$sequential_outcome_regression <- vector(mode = "list",length(run_these_targets))
    names(x$sequential_outcome_regression) = run_these_targets
    for (target_name in run_these_targets){
        message("Running target: ",target_name)
        x$sequential_outcome_regression[[target_name]] <- vector(mode = "list",3)
        names(x$sequential_outcome_regression[[target_name]]) <- c("predicted_values","fit","intervened_data")
        ## protocols <- names(x$protocols)
        protocols <- x$targets[[target_name]]$protocols
        for (protocol_name in protocols){
            message("Current protocol: ",protocol_name)
            protocol_data <- x$prepared_data$data
            N <- NROW(protocol_data)
            # initialize influence curve vector
            for (th in 1:length(time_horizon)){
                x$IC[[target_name]][[protocol_name]][[name_time_horizon[[th]]]] <- numeric(N)
            }
            #
            # g-part: fit nuisance parameter models for propensity and censoring
            #
            current_protocol <- x$protocols[[protocol_name]]
            intervention_table <- current_protocol$intervention_table
            if (length(x$times)>1){
                protocol_Anodes <- sapply(x$times[-length(x$times)],function(tk){
                    paste0(x$protocols[[protocol_name]]$treatment_variables,"_",tk)
                })
            } else{
                protocol_Anodes <- x$protocols[[protocol_name]]$treatment_variables
            }
            protocol_Cnodes <- x$prepared_data$Cnodes
            protocol_Dnodes <- x$prepared_data$Dnodes
            protocol_Ynodes <- x$prepared_data$Ynodes
            #
            # define a matrix with the cumulative intervention/censoring probabilities
            #
            if (refit || length(x$cumulative_intervention_probs[[protocol_name]]) == 0){
                intervention_probs <- data.table(ID = protocol_data[[x$names$id]])
                intervention_probs <- cbind(intervention_probs,matrix(1,nrow = N,ncol = length(protocol_Anodes)+length(protocol_Cnodes)))
                setnames(intervention_probs,c(x$names$id,c(rbind(protocol_Anodes,protocol_Cnodes))))
                # predict the propensity score/1-probability of censored
                # intervene according to protocol for targets
                # in the last time interval we do not need propensities/censoring probabilities  
                for (j in x$times[-length(x$times)]){
                    if (j > 0){
                        outcome_free <- protocol_data[[protocol_Ynodes[[j]]]]%in%0
                        # competing risks
                        if (length(protocol_Dnodes)>0) outcome_free <- outcome_free&protocol_data[[protocol_Dnodes[[j]]]]%in%0
                        uncensored <- protocol_data[[protocol_Cnodes[[j]]]]%in%"uncensored"
                    }else{
                        # assume that all are at risk at time 0
                        outcome_free <- rep(TRUE,N)
                        uncensored <- rep(TRUE,N)
                    }
                    if (any(outcome_free & uncensored)){
                        # fit the propensity regression model
                        if (refit || length(x$models[[protocol_name]][["propensity"]][[protocol_Anodes[[j+1]]]]$fit) == 0){
                            if (!is.null(ff <- x$models[[protocol_name]][["propensity"]][[protocol_Anodes[[j+1]]]]$formula)){
                                x$models[[protocol_name]][["propensity"]][[protocol_Anodes[[j+1]]]]$fit <- do.call(learner,list(formula = ff,data = protocol_data[outcome_free&uncensored],...))
                            }
                        }
                        # fit censoring model
                        if (refit || length(x$models[[protocol_name]][["censoring"]]$fit) == 0){
                            if (!is.null(ff <- x$models[[protocol_name]][["censoring"]][[protocol_Cnodes[[j+1]]]]$formula)){
                                x$models[[protocol_name]][["censoring"]][[protocol_Cnodes[[j+1]]]]$fit <- do.call(learner,list(formula = ff,data = protocol_data[outcome_free&uncensored],...))
                            }
                        }
                    }
                    # set the treatment variables to their protocol values
                    ## FIXME: could subset data to the variables in the current formula
                    # FIXME: also those censored?
                    intervened_data = intervene(data = protocol_data[outcome_free],
                                                intervention_table = intervention_table,
                                                time = j)
                    intervention_probs[outcome_free][[protocol_Anodes[[j+1]]]] <- riskRegression::predictRisk(x$models[[protocol_name]][["propensity"]][[protocol_Anodes[[j+1]]]]$fit,
                                                                                                              newdata = intervened_data)
                    intervention_probs[outcome_free][[protocol_Cnodes[[j+1]]]] <- riskRegression::predictRisk(x$models[[protocol_name]][["censoring"]][[protocol_Cnodes[[j+1]]]]$fit,
                                                                                                              newdata = intervened_data)
                }
                # FIXME: write this rowCumprods in armadillo
                x$cumulative_intervention_probs[[protocol_name]] <- matrixStats::rowCumprods(as.matrix(intervention_probs[,-1,with = FALSE]))
            }
            #
            # define a matrix which indicates if the intervention is followed
            #
            if (length(intervention_match <- x$intervention_match[[protocol_name]]) == 0){
                # fixme for target in targets
                intervention_match <- matrix(0,ncol = length(protocol_Anodes),nrow = N)
                for(j in x$times[-c(length(x$times))]){
                    if (j == 0)
                        intervention_match[,j+1] <- previous <- (protocol_data[[intervention_table[j+1]$variable]] %in% c(intervention_table[j+1]$value,NA))
                    else
                        intervention_match[,j+1] <- previous <- previous*(protocol_data[[intervention_table[j+1]$variable]] %in% c(intervention_table[j+1]$value,NA))
                }
                colnames(intervention_match) <- protocol_Anodes
                x$intervention_match[[protocol_name]] <- intervention_match
            }
            # 
            # Q-part: loop backwards in time through iterative condtional expectations
            #
            x$estimate[[target_name]][[protocol_name]] <- data.table(Target = target_name,
                                                                     Protocol = protocol_name,
                                                                     "Time_horizon" = time_horizon,
                                                                     Estimator = x$targets[[target_name]]$estimator,
                                                                     Estimate = numeric(length(time_horizon)))
            # loop across time-horizons
            for (th in time_horizon){
                x <- sequential_regression(x = x,
                                           target_name = target_name,
                                           protocol_name = protocol_name,
                                           time_horizon = th,
                                           learner = learner,...)
            }
        }
    }
    x
}
######################################################################
### run_rtmle.R ends here
