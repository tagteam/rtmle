### intervention_probabilities.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 17 2024 (09:26) 
## Version: 
## Last-Updated: Oct 28 2024 (16:02) 
##           By: Thomas Alexander Gerds
##     Update #: 32
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
            # see who is at_risk at time_j
            outcome_free_and_uncensored <- x$followup$last_interval >= j
            if (any(outcome_free_and_uncensored)){
                # set the treatment variables to their protocol values
                ## FIXME: could subset data to the variables in the current formula
                ## FIXME: predict also subjects who will be censored in current interval? Think: YES
                intervened_data <- intervene(data = x$prepared_data[outcome_free_and_uncensored][,(1:(-1+match(treatment_variables[[j+1]],names(x$prepared_data)))),with = FALSE],
                                            intervention_table = intervention_table,
                                            time = j)
                for (G in c(treatment_variables[[j+1]],censoring_variables[[j+1]]))
                    # fit the propensity and censoring regression models
                    # and store probabilities as intervention_probs
                    if (refit || length(x$models[[protocol_name]][[G]]$fit) == 0){
                        if (!is.null(ff <- x$models[[protocol_name]][[G]]$formula)){
                            predicted_values <- do.call(learner,list(character_formula = ff,data = x$prepared_data[outcome_free_and_uncensored],intervened_data = intervened_data))
                            x$model[[protocol_name]][[G]]$fit <- attr(predicted_values,"fit")
                            data.table::setattr(predicted_values,"fit",NULL)
                            intervention_probs[outcome_free_and_uncensored][[G]] <- predicted_values
                            # FIXME: this hack works but only when there are exactly 2 treatment levels!
                            if (!(G %in% censoring_variables)){ # then G is a treatment variable
                                if (match(intervention_table[time == j,value],x$names$treatment_levels)<length(x$names$treatment_levels))
                                    intervention_probs[outcome_free_and_uncensored][[G]] <- 1-intervention_probs[outcome_free_and_uncensored][[G]]
                            }
                        }
                    }
            }
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
