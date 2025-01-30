### intervention_probabilities.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Oct 17 2024 (09:26)
## Version:
## Last-Updated: Nov  25 2024 (08:58)
##           By: Alessandra
##     Update #: 100
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
  # FIXME: consider using x$sample_size?
  N <- NROW(x$prepared_data)
  #
  # a matrix with the cumulative intervention/censoring probabilities
  #
  current_protocol <- x$protocols[[protocol_name]]
  intervention_table <- current_protocol$intervention_table
  # added this to check which type of intervention we are doing
  intervention_type <- current_protocol$intervention_type

  # since time_horizon can be a vector we need the maximum
  if (missing(time_horizon))
    max_time_horizon <- max(x$time)
  else
    max_time_horizon <- max(time_horizon)
  eval_times <- x$times[-length(x$times)]
  ## if max_time_horizon = 5 then we need propensities up to time 4
  eval_times <- eval_times[eval_times < max_time_horizon]
  if (length(x$times)>1){
    treatment_variables <- sapply(eval_times,function(tk){
      paste0(x$protocols[[protocol_name]]$treatment_variables,"_",tk)
    })
  } else{
    treatment_variables <- x$protocols[[protocol_name]]$treatment_variables
  }
  censoring_variables <- paste0(x$names$censoring,"_",1:max_time_horizon)
  competing_variables <- paste0(x$names$competing,"_",1:(max_time_horizon-1))
  outcome_variables <- paste0(x$names$outcome,"_",1:max_time_horizon)

  if(intervention_type=="stochastic"){
    ### MODIFY THIS TO HAVE ONLY THE TREATMENT VARIABLE
    # for less memory
    # also, we might have to define this inside the refit
    stochastic_probs <- data.table(ID = x$prepared_data[[x$names$id]])
    stochastic_probs <- cbind(stochastic_probs,matrix(1,nrow = N,ncol = length(treatment_variables)+length(censoring_variables) ))
    setnames(stochastic_probs,c(x$names$id,c(rbind(treatment_variables,censoring_variables))))
  }
  if (refit || length(x$cumulative_intervention_probs[[protocol_name]]) == 0){
    intervention_probs <- data.table(ID = x$prepared_data[[x$names$id]])
    intervention_probs <- cbind(intervention_probs,matrix(1,nrow = N,ncol = length(treatment_variables)+length(censoring_variables) )) # initialize this to a matrix of 1s
    setnames(intervention_probs,c(x$names$id,c(rbind(treatment_variables,censoring_variables))))
    # predict the propensity score/1-probability of censored
    # intervene according to protocol for targets
    # in the last time interval we do not need propensities/censoring probabilities
    # we can initialize to be equal to 1 because we take track of the censored info at each time and we do the average onlu over the uncensored and outcome free
    for (j in eval_times){
      # see who is at_risk at time_j
      outcome_free_and_uncensored <- x$followup$last_interval >= j
      if (any(outcome_free_and_uncensored)){
        if (length(censoring_variables[[j+1]])>0)
          history_of_variables <- 1:(match(censoring_variables[[j+1]],names(x$prepared_data)))
        else
          history_of_variables <- 1:(match(treatment_variables[[j+1]],names(x$prepared_data)))
        # FIXME: would be better to restrict to the variables that occur in the current formula
        current_data <- x$prepared_data[outcome_free_and_uncensored,history_of_variables,with = FALSE]
        # in current data we only have people AT RISK! so we calculate the propensity score only on these (obviously)
        current_constants <- sapply(current_data, function(x){length(unique(x))==1})
        if (any(current_constants)) {
          current_constants <- names(current_constants[current_constants])
        }else{
          current_constants <- NULL
        }
        # set the treatment variables to their protocolled values
        # for the stochastic intervention we still need that for the nuisance parameter of the propensity score
        # in fact we need the predicted value at this time based on the history (depending on the formula)
        # this is why we ll let value as their observed value through intervene_function
        # we should then use the intervene
        if(intervention_type=="stochastic"){
          # intervened data are basically the current_data for the prediction in the stochastic intervention
          intervened_data <- current_data
        }
        else{
          # WHY do we need this? cause we do not need to predict the propensity score in the intervened data, not even for the statics intervention
          # also at the end you ll calculate the average only on the ones that followed the intervention, thx to intervention_match table
          intervened_data <- do.call(x$protocol[[protocol_name]]$intervene_function,
                                     list(data = current_data,
                                          intervention_table = intervention_table,
                                          time = j)) # not sure why time is there, but maybe will make sense
        }


        for (G in c(treatment_variables[[j+1]],censoring_variables[[j+1]]))
          # fit the propensity and censoring regression models
          # and store probabilities as intervention_probs
          if (refit || length(x$models[[protocol_name]][[G]]$fit) == 0){
            if (G %in% current_constants){
              if (G %in% censoring_variables){
                predicted_values <- 1*(current_data[[G]] == x$names$uncensored_label)
              }else{
                predicted_values <- current_data[[G]] ### WHY THIS?
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

                # in args we define the argument for the learner
                # where we have intervened data which are the current data with the value for the treatment needed for that specific intervention
                # for the stochastic one it is basically the current_data
                args <- list(character_formula = ff,
                             data = current_data[,!(names(current_data)%in%current_constants),with = FALSE],
                             intervened_data = intervened_data[,!(names(intervened_data) %in% current_constants),with = FALSE])

                if (length(learner)>1){
                  args <- c(args,
                            list(learners = learner,
                                 outcome_variable = G,
                                 id_variable = x$names$id))
                  if (inherits(try(
                    predicted_values <- do.call("superlearn",c(args,list(seed = seed))),silent = FALSE),
                    "try-error")) {
                    ## browser(skipCalls=1L)
                    stop(paste0("Failed to superlearn/crossfit with formula ",ff))
                  }
                }else{
                  if (inherits(try(
                    predicted_values <- do.call(learner,args),silent = FALSE),
                    "try-error")) {
                    ## browser(skipCalls=1L)
                    stop(paste0("Failed to learn/predict with formula ",ff))
                  }
                }
              }
              x$models[[protocol_name]][[G]]$fit <- attr(predicted_values,"fit")
              data.table::setattr(predicted_values,"fit",NULL)
              intervention_probs[outcome_free_and_uncensored][[G]] <- predicted_values # this is where I have \hat{g_t} (propensity score saved)

              # we have to calculate the stochastic probability in the observed data (current_data)
              if(intervention_type=="stochastic" & G %in% treatment_variables){
              stochastic_probs[outcome_free_and_uncensored][[G]]<-do.call(x$protocol[[protocol_name]]$intervene_function,
                                                                          list(data = current_data,
                                                                               intervention_table=intervention_table,
                                                                               current.time = j))
              }
              # FIXME: this hack works but only when there are exactly 2 treatment levels! (what is this for?)
              if (!(G %in% censoring_variables)){ # then G is a treatment variable
                #  I do not see how this would be needed in the stochastic intervention
                # however, we have to calculate the ratio as the weights:
                # to do the cumulative product afterwards

                if(intervention_type!="stochastic"){
                # NOT SURE WHAT THIS IS DOING, for a stochastic intervention we do not have to do this step because we always need \hat{p}
                if (match(intervention_table[time == j,value],x$names$treatment_levels)<length(x$names$treatment_levels))
                  intervention_probs[outcome_free_and_uncensored][[G]] <- 1-intervention_probs[outcome_free_and_uncensored][[G]]
                }
              }
            }
          }
      }
    }
    # FIXME: remove this when not needed anymore
    x$intervention_probs[[protocol_name]] <- intervention_probs # here we have th propensity score for treatment and censoring
    # FIXME: write this rowCumprods in armadillo
    x$cumulative_intervention_probs[[protocol_name]] <- matrixStats::rowCumprods(as.matrix(intervention_probs[,-1,with = FALSE]))

  }
  #
  # define a matrix which indicates if the intervention is followed
  # in the stochastic intervention there is not such a thing,
  # we have to calculate the stochastic intervention in the observed data instead
  # and we can do that before, no need for another loop
  #
  if(intervention_type!="stochastic"){
  if (length(intervention_match <- x$intervention_match[[protocol_name]]) == 0
      ||
      ## the previous run could have produced the matrix but maybe not for all time points
      NCOL(x$intervention_match[[protocol_name]])<length(eval_times)){
    intervention_match <- matrix(0,ncol = length(treatment_variables),nrow = N) # this only refers to the treatment "match" to the target intervention
    for(j in eval_times){
      if (j == 0)
        intervention_match[,j+1] <- previous <- (x$prepared_data[[intervention_table[j+1]$variable]] %in% c(intervention_table[j+1]$value,NA))
      else
        intervention_match[,j+1] <- previous <- previous*(x$prepared_data[[intervention_table[j+1]$variable]] %in% c(intervention_table[j+1]$value,NA))
    }
    colnames(intervention_match) <- treatment_variables
    x$intervention_match[[protocol_name]] <- intervention_match
  }
  }
  else{

    x$stochastic_probs[[protocol_name]] <- stochastic_probs # here we have the estimated stochastic probability in the observed data
    x$cumulative_stochastic_probs[[protocol_name]] <- matrixStats::rowCumprods(as.matrix(stochastic_probs[,-1,with = FALSE]))
    # here we have to divide the stochastic intervention with the propensity score of treatment and censoring:
    # stochastic weights (here we have basically the definition of the clever covariate)
    x$cumulative_intervention_probs[[protocol_name]]<- x$cumulative_stochastic_probs[[protocol_name]]/ x$cumulative_intervention_probs[[protocol_name]]

    intervention_match <- matrix(1,ncol = length(treatment_variables),nrow = N)
    colnames(intervention_match) <- treatment_variables
    x$intervention_match[[protocol_name]] <-intervention_match # we put 1 so that in the sequential regression considers everyone at risk for the weight (all matches)
  }
  x
}

######################################################################
### intervention_probabilities.R ends here
