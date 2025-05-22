### intervention_probabilities.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Oct 17 2024 (09:26)
## Version:
## Last-Updated: May 21 2025 (15:49)
##           By: Alessandra
##     Update #: 221
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
  ##
  treatment_variables <- lapply(eval_times,function(this_time){
    x$protocols[[protocol_name]]$intervention_table[time == this_time]$variable
  })
  if (length(x$names$censoring)>0){
    censoring_variables <- paste0(x$names$censoring,"_",1:max_time_horizon)
  }else{
    censoring_variables <- NULL
  }

  ## if there is competing event
  if (length(x$names$competing)>0){
    competing_variables <- paste0(x$names$competing,"_",1:(max_time_horizon-1))
  }else{
    competing_variables <- NULL
  }


  outcome_variables <- paste0(x$names$outcome,"_",1:max_time_horizon)

  # order vector of treatment and censoring_variables
  G_names <- unlist(lapply(eval_times,function(this_time){
    c(treatment_variables[[this_time+1]],censoring_variables[[this_time+1]])
  }))
  if (refit || NCOL(x$cumulative_intervention_probs[[protocol_name]]) < length(G_names)){


    ## when we define a stochastic intervention we have to calculate the
    # stochastic probability of the treatment at each time
    if(intervention_type=="stochastic"){
      ### MODIFY THIS TO HAVE ONLY THE TREATMENT VARIABLE
      # for less memory
      stochastic_probs <- data.table(ID = x$prepared_data[[x$names$id]])
      stochastic_probs <- cbind(stochastic_probs,matrix(1,
                                                        nrow = N,
                                                        ncol =  length(unlist(treatment_variables))+length(censoring_variables) ))
      setnames(stochastic_probs,new = c(x$names$id,G_names))
    }


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
        current_G_variables <- c(treatment_variables[[j+1]],censoring_variables[[j+1]])
        ## for dynamic or stochastic intervention, I guess that all the needed variables
        ## will be included in the propsensity score model, so we do not need to add other stuff!
        ## if this is not the case we should understand how to get them!
        current_formulas <- lapply(current_G_variables,function(G){
          formula(x$models[[protocol_name]][[G]]$formula)
        })
        names(current_formulas) <- current_G_variables
        # all variables that are either outcome or predictor in this time interval
        history_of_variables <- unique(unlist(sapply(current_formulas,all.vars)))
        # data used to fit the models in this time interval
        current_data <- x$prepared_data[outcome_free_and_uncensored,history_of_variables,with = FALSE]
        # checking for missing values in any of the predictor variables
        if (any(is.na(current_data))){
          has_missing <- sapply(current_data,function(x)sum(is.na(x)))
          has_missing <- has_missing[has_missing != 0]
          stop(paste0("Missing values detected in data for fitting nuisance parameter models at time ",j,":\n",
                      paste0(names(has_missing),": n=",has_missing,collapse = "\n")))
        }
        current_constants <- sapply(current_data, function(x){length(unique(x))==1})
        if (any(current_constants)) {
          current_constants <- names(current_constants[current_constants])
        }else{
          current_constants <- NULL
        }

        if(intervention_type=="stochastic"){
          # intervened data are basically the current_data (observed)
          #  for the prediction in the stochastic intervention
          intervened_data <- current_data
        }
        else{
        # set the treatment variables to their protocolled values
        # this step is useful for the dynamic intervention later (as the function is constructed so far!)
        intervened_data <- do.call(x$protocol[[protocol_name]]$intervene_function,
                                   list(data = current_data,
                                        intervention_table = intervention_table,
                                        current.time = j))
        }
        for (G in current_G_variables){
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

              # FIXME: this hack works but only when there are exactly 2 treatment levels!
              if (!(G %in% censoring_variables)){ # then G is a treatment variable
                levels <- x$names$treatment_levels
                if(intervention_type=="stochastic"){
                  # if stochastic intervention we have the deinition in intervene_function with the value (P(A=1| X)) where o is the ref level for the factor A
                  # so for people that were not treated, we have to give 1-prob
                  # and we have to change it for both, propensity score(predicted_values) and stochastic probs
                  # based on the observed data:
                  prob_stochastics<-do.call(x$protocol[[protocol_name]]$intervene_function,
                                            list(data = current_data,
                                                 name.treat=G, ## the treatment we want to calculate the prob for (because in case of double treatment, only time=j will not be sufficient)
                                                 current.time = j))


                  if(length(levels)==2){
                    ## in case of categorical treatment, the function will provide a matrix with prob for
                    ## the obsevred avlue in the current data, we could actually ask to do it also for the binary case:
                    ## FIX ME: this would work only with 0 or 1
                  ind_trt0<-which(intervened_data[, ..G ]==levels[1]) # ref level in the specific treatment
                  prob_stochastics[ind_trt0]<-1-prob_stochastics[ind_trt0]
                  predicted_values[ind_trt0] <- 1-predicted_values[ind_trt0]
                  stochastic_probs[outcome_free_and_uncensored][[G]]<-prob_stochastics
                  }

                }

                else{
                  # this would work for both the dynamic intervention and the static one
                  # in the dynamic one the intervened_data are defined by the intervene_function
                  # that provides the treatment value based on some rules.

                  if(length(x$names$treatment_levels)>2){
                    predicted_values<-predicted_values[,match(intervention_table[time == j,value],x$names$treatment_levels)]}
                  else{
                  ind_trt0<-(intervened_data[, ..G] == x$names$treatment_levels[1])
                  predicted_values[ind_trt0] <- 1 - predicted_values[ind_trt0]
                  }

                # FIXME: could extract levels from data or (even better) add an association list
                #        list(A=c("A_0","A_1","A_2"), B=c("B_0","B_1","B_2"),... to the object x
                # if (length(levels) == 2){
                #   if (match(intervention_table[variable == G][["value"]],levels) == 1) ##i s the reference level
                #     # probs are 1-probs
                #     intervention_probs[outcome_free_and_uncensored][[G]] <- 1-intervention_probs[outcome_free_and_uncensored][[G]]
                # }

            }
          }

          # I moved this down cause predicted_values was a matrix if for a multinomial, and we only need one
          intervention_probs[outcome_free_and_uncensored][[G]] <- predicted_values
           }
         }
        } # close the for over G at the current time
        } # close the if(any at risk)
     } # close the for over time
    # FIXME: remove this when not needed anymore
    x$intervention_probs[[protocol_name]] <- intervention_probs
    # FIXME: write this rowCumprods in armadillo
    x$cumulative_intervention_probs[[protocol_name]] <- matrixStats::rowCumprods(as.matrix(intervention_probs[,-1,with = FALSE]))
  }
  #


  if(intervention_type!="stochastic"){

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

    if(intervention_type=="dynamic"){
      ## for the dynamic intervention we have to check the observed value with the intervened value:
      ## in case of multiple treatment, then this will intervene in both of them
      ## we can also construct this in the previous loop
      intervention_values_tab <- do.call(x$protocol[[protocol_name]]$intervene_function,
                                         list(data = x$prepared_data,
                                              intervention_table = intervention_table))

      for(j in eval_times){
        intervention_variables <- intervention_table[time == j]$variable
        observed_values <- x$prepared_data[,intervention_variables,with = FALSE]
        for (v in 1:length(intervention_variables)){ # this for is needed in case of multiple intervention variables
          intervention_match[,j+1] <- previous <- previous*(intervention_values_tab[[intervention_variables[[v]]]] - observed_values[[intervention_variables[[v]]]] ==0)
        }

      }
    }

    else{
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
    }
    ## the intervention_match matrix has one column per time point
    # when there are mulitple treatment variables we paste-collapse the names
    colnames(intervention_match) <- sapply(treatment_variables,paste,collapse = ",")
    x$intervention_match[[protocol_name]] <- intervention_match
   }
  }
      else{
        x$stochastic_probs[[protocol_name]] <- stochastic_probs # here we have the estimated stochastic probability in the observed data
        x$cumulative_stochastic_probs[[protocol_name]] <- matrixStats::rowCumprods(as.matrix(stochastic_probs[,-1,with = FALSE]))
        # here we have to divide the stochastic intervention with the propensity score of treatment and censoring:
        # stochastic weights (here we have basically the definition of the clever covariate)
        x$cumulative_intervention_probs[[protocol_name]]<- x$cumulative_stochastic_probs[[protocol_name]]/ x$cumulative_intervention_probs[[protocol_name]]
        intervention_match <- matrix(1,ncol = length(eval_times),nrow = N)
        # when there are multiple treatment variables we paste-collapse the names
        colnames(intervention_match) <- sapply(treatment_variables,paste,collapse = ",")
        x$intervention_match[[protocol_name]] <-intervention_match # we put 1 so that in the sequential regression considers everyone at risk for the weight (all matches)
      }


  x
}

######################################################################
### intervention_probabilities.R ends here
