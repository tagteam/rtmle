### intervention_probabilities.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 17 2024 (09:26) 
## Version: 
## Last-Updated: Oct 19 2024 (09:39) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
intervention_probabilities <- function(x,
                                       protocol_name,
                                       refit = FALSE,
                                       learner,
                                       ...){
    # FIXME: consider using x$sample_size?
    N <- NROW(x$prepared_data)
    #
    # a matrix with the cumulative intervention/censoring probabilities
    #
    current_protocol <- x$protocols[[protocol_name]]
    intervention_table <- current_protocol$intervention_table
    if (length(x$times)>1){
        treatment_variables <- sapply(x$times[-length(x$times)],function(tk){
            paste0(x$protocols[[protocol_name]]$treatment_variables,"_",tk)
        })
    } else{
        treatment_variables <- x$protocols[[protocol_name]]$treatment_variables
    }
    # since time_horizon can be a vector we need the maximum
    max_time_horizon <- max(x$time)
    censoring_variables <- paste0(x$names$censoring,"_",1:max_time_horizon)
    competing_variables <- paste0(x$names$competing,"_",1:(max_time_horizon-1))
    outcome_variables <- paste0(x$names$outcome,"_",1:max_time_horizon)
    if (refit || length(x$cumulative_intervention_probs[[protocol_name]]) == 0){
        intervention_probs <- data.table(ID = x$prepared_data[[x$names$id]])
        intervention_probs <- cbind(intervention_probs,matrix(1,nrow = N,ncol = length(treatment_variables)+length(censoring_variables)))
        setnames(intervention_probs,c(x$names$id,c(rbind(treatment_variables,censoring_variables))))
        # predict the propensity score/1-probability of censored
        # intervene according to protocol for targets
        # in the last time interval we do not need propensities/censoring probabilities
        for (j in x$times[-length(x$times)]){
            if (j > 0){
                outcome_free <- x$prepared_data[[outcome_variables[[j]]]]%in%0
                # competing risks
                if (length(competing_variables)>0) outcome_free <- outcome_free&x$prepared_data[[competing_variables[[j]]]]%in%0
                uncensored <- x$prepared_data[[censoring_variables[[j]]]]%in%x$names$uncensored_label
            }else{
                # assume that all are at risk at time 0
                outcome_free <- rep(TRUE,N)
                uncensored <- rep(TRUE,N)
            }
            if (any(outcome_free & uncensored)){
                # fit the propensity regression model
                if (
                    refit || length(x$models[[protocol_name]][["propensity"]][[treatment_variables[[j+1]]]]$fit) == 0){
                    if (!is.null(ff <- x$models[[protocol_name]][["propensity"]][[treatment_variables[[j+1]]]]$formula)){
                        x$models[[protocol_name]][["propensity"]][[treatment_variables[[j+1]]]]$fit <- do.call(learner,list(formula = ff,data = x$prepared_data[outcome_free&uncensored],...))
                    }
                }
                # fit censoring model
                if (refit || length(x$models[[protocol_name]][["censoring"]]$fit) == 0){
                    if (!is.null(ff <- x$models[[protocol_name]][["censoring"]][[censoring_variables[[j+1]]]]$formula)){
                        x$models[[protocol_name]][["censoring"]][[censoring_variables[[j+1]]]]$fit <- do.call(learner,list(formula = ff,data = x$prepared_data[outcome_free&uncensored],...))
                    }
                }
            }
            # set the treatment variables to their protocol values
            ## FIXME: could subset data to the variables in the current formula
            ## FIXME: predict also subjects who will be censored in current interval? Think: YES
            intervened_data = intervene(data = x$prepared_data[outcome_free],
                                        intervention_table = intervention_table,
                                        time = j)
            intervention_probs[outcome_free][[treatment_variables[[j+1]]]] <- riskRegression::predictRisk(x$models[[protocol_name]][["propensity"]][[treatment_variables[[j+1]]]]$fit,
                                                                                                          newdata = intervened_data)
            # FIXME: this hack works probably only with exactly 2 treatment levels?
            if (match(intervention_table[time == j,value],x$names$treatment_levels)<length(x$names$treatment_levels))
                intervention_probs[outcome_free][[treatment_variables[[j+1]]]] <- 1-intervention_probs[outcome_free][[treatment_variables[[j+1]]]]
            intervention_probs[outcome_free][[censoring_variables[[j+1]]]] <- riskRegression::predictRisk(x$models[[protocol_name]][["censoring"]][[censoring_variables[[j+1]]]]$fit,
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
        intervention_match <- matrix(0,ncol = length(treatment_variables),nrow = N)
        for(j in x$times[-c(length(x$times))]){
            if (j == 0)
                intervention_match[,j+1] <- previous <- (x$prepared_data[[intervention_table[j+1]$variable]] %in% c(intervention_table[j+1]$value,NA))
            else
                intervention_match[,j+1] <- previous <- previous*(x$prepared_data[[intervention_table[j+1]$variable]] %in% c(intervention_table[j+1]$value,NA))
        }
        colnames(intervention_match) <- treatment_variables
        x$intervention_match[[protocol_name]] <- intervention_match
    }
    x
}

######################################################################
### intervention_probabilities.R ends here
