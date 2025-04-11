### intervention_probabilities.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 17 2024 (09:26) 
## Version: 
## Last-Updated: Apr 11 2025 (16:39) 
##           By: Thomas Alexander Gerds
##     Update #: 206
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
                                       time_horizon,
                                       seed,
                                       ...){
    variable <- NULL
    # FIXME: consider using x$sample_size?
    N <- NROW(x$prepared_data)
    #
    # a matrix with the cumulative intervention/censoring probabilities
    #
    current_protocol <- x$protocols[[protocol_name]]
    intervention_table <- current_protocol$intervention_table
    # since time_horizon can be a vector we need the maximum
    if (missing(time_horizon))
        max_time_horizon <- max(x$time)
    else
        max_time_horizon <- max(time_horizon)
    eval_times <- x$times[-length(x$times)]
    ## if max_time_horizon = 5 then we need propensities up to time 4
    eval_times <- eval_times[eval_times < max_time_horizon]
    ##
    treatment_variables <- lapply(eval_times,function(this_time){
        x$protocols[[protocol_name]]$intervention_table[time == this_time]$variable
    })
    if (length(x$names$censoring)>0){
        censoring_variables <- paste0(x$names$censoring,"_",1:max_time_horizon)
    }else{
        censoring_variables <- NULL
    }
    competing_variables <- paste0(x$names$competing,"_",1:(max_time_horizon-1))
    outcome_variables <- paste0(x$names$outcome,"_",1:max_time_horizon)
    # order vector of treatment and censoring_variables
    G_names <- unlist(lapply(eval_times,function(this_time){
        c(treatment_variables[[this_time+1]],censoring_variables[[this_time+1]])
    }))
    if (refit || NCOL(x$cumulative_intervention_probs[[protocol_name]]) < length(G_names)){
        # FIXME: the following code should be improved to increase the readability
        intervention_probs <- data.table(ID = x$prepared_data[[x$names$id]])
        intervention_probs <- cbind(intervention_probs,
                                    matrix(1,
                                           nrow = N,
                                           ncol = length(unlist(treatment_variables))+length(censoring_variables)))
        setnames(intervention_probs,new = c(x$names$id,G_names))
        # predict the propensity score/1-probability of censored
        # intervene according to protocol for targets
        # in the last time interval we do not need propensities/censoring probabilities
        for (j in eval_times){
            # see who is at_risk at time_j
            if (j == 0){
                outcome_free_and_uncensored <- rep(TRUE,N)
            }else{
                outcome_free_and_uncensored <- x$followup$last_interval >= j
            }
            if (any(outcome_free_and_uncensored)){
                if (length(censoring_variables[[j+1]])>0){
                    last_G <- censoring_variables[[j+1]]
                } else{
                    # FIXME: this is not easy to read, but when there are multiple treatment variables
                    #        at time j, such as A_j and B_j then we match on the last
                    last_G <- rev(treatment_variables[[j+1]])[[1]]
                }
                # FIXME: this requires that the data are chronologically ordered
                history_of_variables <- 1:(match(last_G,names(x$prepared_data)))
                # FIXME: would be better to restrict to the variables that occur in the current formula
                # FIXME: would be better to remove the id variable
                current_data <- x$prepared_data[outcome_free_and_uncensored,history_of_variables,with = FALSE]
                current_constants <- sapply(current_data, function(x){length(unique(x))==1})
                if (any(current_constants)) {
                    current_constants <- names(current_constants[current_constants])
                }else{
                    current_constants <- NULL
                }
                # set the treatment variables to their protocolled values
                intervened_data <- do.call(x$protocol[[protocol_name]]$intervene_function,
                                           list(data = current_data,
                                                intervention_table = intervention_table,
                                                time = j))
                for (G in c(treatment_variables[[j+1]],censoring_variables[[j+1]])){
                    # fit the propensity and censoring regression models
                    # and store probabilities as intervention_probs
                    if (refit || length(x$models[[protocol_name]][[G]]$fit) == 0){
                        if (G %in% current_constants){
                            if (G %in% censoring_variables){
                                predicted_values <- 1*(current_data[[G]] == x$names$uncensored_label)
                            }else{
                                predicted_values <- current_data[[G]]
                            }
                            x$models[[protocol_name]][[G]]$fit <- "No variation in this variable"
                        }else{
                            if (length(ff <- x$models[[protocol_name]][[G]]$formula)>0){
                                # remove constant predictor variables
                                ff_vars <- all.vars(stats::formula(ff))
                                if (length(current_constants)>0){
                                    ## warning("Removing constant variables at time ",j,":\n",paste0(current_constants,collapse = ", ") )
                                    ff <- delete_variables_from_formula(character_formula = ff,delete_vars = current_constants)
                                }
                                args <- list(character_formula = ff,
                                             data = current_data[,!(names(current_data)%in%current_constants),with = FALSE],
                                             intervened_data = intervened_data,...)
                                if (length(learner)>1){
                                    args <- c(args,
                                              list(learners = learner,
                                                   outcome_variable = G,
                                                   outcome_target_level = levels(current_data[[G]])[[2]],
                                                   id_variable = x$names$id))
                                    if (inherits(try(
                                        predicted_values <- do.call("superlearn",c(args,list(seed = seed))),silent = FALSE),
                                        "try-error")) {
                                        stop(paste0("Failed to superlearn/crossfit with formula ",ff))
                                    }
                                }else{
                                    if (inherits(try(
                                        predicted_values <- do.call(learner,args),silent = FALSE),
                                        "try-error")) {
                                        stop(paste0("Failed to learn/predict with formula ",ff))
                                    }
                                }
                            }else{
                                stop(paste0("Cannot see a formula for estimating regression of ",G," at x$models[['",protocol_name,"']][['",G,"']]$formula"))
                            }
                            x$models[[protocol_name]][[G]]$fit <- attr(predicted_values,"fit")
                            data.table::setattr(predicted_values,"fit",NULL)
                            intervention_probs[outcome_free_and_uncensored][[G]] <- predicted_values
                            # FIXME: this hack works but only when there are exactly 2 treatment levels!
                            ## if (G == "B_1") browser(skipCalls=1L)
                            ## if (G == "A_1") browser(skipCalls=1L)
                            if (!(G %in% censoring_variables)){ # then G is a treatment variable
                                # FIXME: could extract levels from data or (even better) add an association list
                                #        list(A=c("A_0","A_1","A_2"), B=c("B_0","B_1","B_2"),... to the object x
                                if (length(levels <- x$names$treatment_options[[sub("_[0-9]+$","",G)]]) == 2){
                                    if (match(intervention_table[variable == G][["value"]],levels) == 1)
                                        # probs are 1-probs
                                        intervention_probs[outcome_free_and_uncensored][[G]] <- 1-intervention_probs[outcome_free_and_uncensored][[G]]
                                }
                            }
                        }
                    }
                }
            }
        }
        # FIXME: remove this when not needed anymore
        x$intervention_probs[[protocol_name]] <- intervention_probs
        # FIXME: write this rowCumprods in armadillo
        x$cumulative_intervention_probs[[protocol_name]] <- matrixStats::rowCumprods(as.matrix(intervention_probs[,-1,with = FALSE]))
    }
    #
    # define a matrix which indicates if the intervention is followed
    # the matrix should have a row for each individual and a column for
    # each time point 
    #
    if (length(intervention_match <- x$intervention_match[[protocol_name]]) == 0
        ||
        ## the previous run could have produced the matrix but maybe not for all time points
        NCOL(x$intervention_match[[protocol_name]])<length(eval_times)){
        intervention_match <- matrix(0,ncol = length(eval_times),nrow = N)
        previous <- rep(1,N)
        for(j in eval_times){
            intervention_variables <- intervention_table[time == j]$variable
            observed_values <- x$prepared_data[,intervention_variables,with = FALSE]
            for (v in 1:length(intervention_variables)){
                # when there are multiple intervention variables
                # all observed values must match  
                intervention_values <- intervention_table[time == j & variable == intervention_variables[[v]]]$value
                intervention_match[,j+1] <- previous <- previous*(observed_values[[intervention_variables[[v]]]] %in% intervention_values)
            }
        }
        ## the intervention_match matrix has one column per time point
        # when there are mulitple treatment variables we paste-collapse the names
        colnames(intervention_match) <- sapply(treatment_variables,paste,collapse = ",")
        x$intervention_match[[protocol_name]] <- intervention_match
    }
    x
}

######################################################################
### intervention_probabilities.R ends here
