### sequential_regression.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Sep 30 2024 (14:30)
## Version:
## Last-Updated: 11 Dec 2024 (08:59)
##           By: Alessandra
##     Update #: 155
#----------------------------------------------------------------------
##
### Commentary: We want to adapt this function to stochastic intervention
#               this implies that we have to calculate the predicted outcome for all possible values of A
#               and we have to use the stochastic intervention probabilities for the expextation of the outcome
#               when using the update_Q we have to specify that it is a stochastic intervention, because the weights are
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
sequential_regression <- function(x,
                                  target_name,
                                  protocol_name,
                                  time_horizon,
                                  learner,
                                  seed = seed,
                                  ...){
  time = Time_horizon = Estimate = Standard_error = Lower = Upper = NULL
  # FIXME: inconsistent listing:
  intervention_table <- x$protocols[[protocol_name]]$intervention_table # this would have value=NULL
  intervention_match <- x$intervention_match[[protocol_name]] ## do not need this for the stochastic intervention
  intervention_type <- x$protocols[[protocol_name]]$intervention_type
  # we might initialize it to null and not change it in the intervention_probabilities function where it has been defined
  if (length(x$times)>1){
    treatment_variables <- sapply(x$times[-length(x$times)],function(tk){
      paste0(x$protocols[[protocol_name]]$treatment_variables,"_",tk)
    })
  } else{
    treatment_variables <- x$protocols[[protocol_name]]$treatment_variables
  }
  censoring_variables <- paste0(x$names$censoring,"_",1:time_horizon)
  competing_variables <- paste0(x$names$competing,"_",1:(time_horizon-1))
  outcome_variables <- paste0(x$names$outcome,"_",1:time_horizon)
  # the first step is the outcome regression of the last time interval
  x$sequential_outcome_regression[[target_name]] = vector(3,mode = "list")
  stopifnot(time_horizon>0)
  label_time_horizon <- paste0("time_horizon_",time_horizon)
  reverse_time_scale <- rev(seq(1,time_horizon,1))
  #start from the last time til the first
  for (j in reverse_time_scale){ # time_horizon is only one number?
    # subjects that are outcome-free and uncensored at the BEGINNING of the interval are 'atrisk'
    outcome_free_and_uncensored <- (x$followup$last_interval >= (j-1))
    # Note that at the time_horizon subjects who are censored in last interval have an NA outcome
    # but at any subsequent (lower time) we can predict the outcome also for subjects who are censored
    # during the current interval. So, while their observed outcome is NA their predicted_outcome is available
    # we keep the NA values and use outcome_free_and_uncensored rather than outcome_free_and_still_uncensored
    # because the TMLE update step needs to subset the predicted values to the non-NA outcomes
    if (j == time_horizon) { #if it is the last time, then the outer expectation:
      outcome_name <- outcome_variables[[time_horizon]]
      # FIXME: protocols could share the outcome formula? YES
      interval_outcome_formula = x$models[[protocol_name]][[outcome_variables[[j]]]]$formula
    } else {
      outcome_name <- "rtmle_predicted_outcome" #otherwise we use what estimated at the previous step
      interval_outcome_formula <- stats::update(stats::formula(x$models[[protocol_name]][[outcome_variables[[j]]]]$formula),"rtmle_predicted_outcome~.")
    }
    ## HERE
    ## Y <- x$prepared_data[outcome_free_and_uncensored][[outcome_name]]
    Y <- x$prepared_data[[outcome_name]]
    # intervene according to protocol for targets
    # FIXME: intervene on all variables or only those after
    #        time j? those in current outcome_formula
    history_of_variables <- names(x$prepared_data)[1:(-1+match(outcome_variables[[j]],names(x$prepared_data)))]
    intervenable_history <- setdiff(history_of_variables,c(outcome_variables,censoring_variables,competing_variables))


    if(intervention_type=="stochastic"){
      # for the stochastic intervention we need to integrate on the all values of A
      # so we have to have the prediction for the intervened data setting A=0 and A=1 (all possible values of A)
      # and respect to the stochastic intervention, so we have to calculate d(1| ) and d(0| )
      intervened_data<- list(data_0=do.call(intervene,
                                list(data = x$prepared_data[,intervenable_history,with = FALSE],
                                     intervention_table = intervention_table[,value:=0],
                                     time = j)) ,
                                data_1=do.call(intervene,
                                  list(data = x$prepared_data[,intervenable_history,with = FALSE],
                                       intervention_table = intervention_table[,value:=1],
                                       time = j)) )

    }
    else{

    intervened_data <- do.call(x$protocol[[protocol_name]]$intervene_function,
                               list(data = x$prepared_data[,intervenable_history,with = FALSE],
                                    intervention_table = intervention_table,
                                    time = j))
    }
    # fit outcome regression
    # we always have the censoring variable in time interval _j *before* the outcome variable in same interval
    # hence to analyse Y_j we need C_j = "uncensored" for the modeling of Y_j
    # we thus remove subjects who are at risk at the beginning of the interval
    # but get censored during the interval (in fact, their outcome is NA)
    if (j == time_horizon)
      learn_variables <- c(history_of_variables,outcome_variables[[j]])
    else
      learn_variables <- c(history_of_variables,"rtmle_predicted_outcome")

    current_data <- x$prepared_data[outcome_free_and_uncensored,learn_variables,with = FALSE]
    current_constants <- sapply(current_data, function(x){length(unique(x))==1})
    if (any(current_constants)) {
      current_constants <- names(current_constants[current_constants])
    }else{
      current_constants <- NULL
    }
    # remove constant predictor variables
    interval_outcome_formula_vars <- all.vars(stats::formula(interval_outcome_formula))
    if (length(current_constants)>0){
      # FIXME: table warnings and show number of times per warning
      ## x$warnings <- paste0("Removing constant variables at time ",j,":\n",paste0(current_constants,collapse = ", "))
      interval_outcome_formula <- delete_variables_from_formula(character_formula = interval_outcome_formula,delete_vars = current_constants)
    }

    if(intervention_type=="stochastic"){

   args<-list(character_formula = interval_outcome_formula,
              data = current_data[,!(names(current_data)%in%current_constants),with = FALSE],
              intervened_data = rbind(intervened_data$data_0[outcome_free_and_uncensored,!(names(intervened_data) %in% current_constants),with = FALSE],
                                      intervened_data$data_1[outcome_free_and_uncensored,!(names(intervened_data) %in% current_constants),with = FALSE]) )
      }

else{
    args <- list(character_formula = interval_outcome_formula,
                 data = current_data[,!(names(current_data)%in%current_constants),with = FALSE],
                 intervened_data = intervened_data[outcome_free_and_uncensored,!(names(intervened_data) %in% current_constants),with = FALSE])}

    if (length(learner)>1){
      if (j == time_horizon)
        args <- c(args,list(learners = learner,outcome_variable = outcome_variables[[j]], id_variable = x$names$id))
      else
        args <- c(args,list(learners = learner,outcome_variable = "rtmle_predicted_outcome", id_variable = x$names$id))
      if (inherits(try(
        fit_last <- do.call("superlearn",c(args,list(seed = seed))),silent = FALSE),
        "try-error")) {
        ## browser(skipCalls=1L)
        stop(paste0("Sequential regression fit failed with formula:\n",interval_outcome_formula))
      }
    }else{
      if (inherits(try(
        fit_last <- do.call(learner,args),silent = FALSE),
        "try-error")){
        ## browser(skipCalls=1L)
        stop(paste0("Sequential regression fit failed with formula:\n",interval_outcome_formula))
      }
    }
    # save fitted object
    x$models[[protocol_name]][[outcome_variables[[j]]]]$fit <- attr(fit_last,"fit")
    data.table::setattr(fit_last,"fit",NULL)
    fit_last <- as.numeric(fit_last) # prediction for the logistic regression, still need the update
    # set predicted value as outcome for next regression
    ## FIXME: if (target$estimator == "tmle")
    old_action <- options()$na.action
    on.exit(options(na.action = old_action))
    options(na.action = "na.pass")
    # note that we can still predict those who are censored at C_j but uncensored at C_{j-1}
    ## y <- riskRegression::predictRisk(fit_last,newdata = intervened_data)
    # avoid missing values due to logit
    if (x$targets[[target_name]]$estimator == "tmle"){
      if (any(fit_last[!is.na(fit_last)] <= 0)) fit_last <- pmax(fit_last,0.0001)
      if (any(fit_last[!is.na(fit_last)] > 1)) fit_last <- pmin(fit_last,0.9999)
    }
    ## y <- pmax(pmin(y,0.99999),0.00001)
    ## y <- pmax(pmin(y,0.9999),0.0001)
    # TMLE update step
    if (length(x$targets[[target_name]]$estimator) == 0 || x$targets[[target_name]]$estimator == "tmle"){
      # use only data from subjects who are uncensored in current interval



      if(intervention_type=="stochastic"){

        # construction of clever covariates
        dim_unc<-sum(outcome_free_and_uncensored)
        Wold <- rep(NA,length(Y))
        Wold[outcome_free_and_uncensored] <- lava::logit(fit_last[1:dim_unc])
        #Y and outcome_free_and_uncensored have double length as for fit_last:
        W0 <- update_Q(Y = Y,
                      logitQ = Wold,
                      cum.g = x$cumulative_intervention_probs[[protocol_name]][,censoring_variables[[j]]],
                      uncensored_undeterministic = outcome_free_and_uncensored,
                      intervention.match = x$intervention_match[[protocol_name]][,intervention_table[time == j-1]$variable], # this would be ones anyways for the stochastic intervention
                      intervention_type=intervention_type)

        Wold <- rep(NA,length(Y))
        Wold[outcome_free_and_uncensored] <- lava::logit(fit_last[(dim_unc+1):length(fit_last)])
        W1 <- update_Q(Y = Y,
                       logitQ = Wold,
                       cum.g = x$cumulative_intervention_probs[[protocol_name]][,censoring_variables[[j]]],
                       uncensored_undeterministic = outcome_free_and_uncensored,
                       intervention.match = x$intervention_match[[protocol_name]][,intervention_table[time == j-1]$variable], # this would be ones anyways for the stochastic intervention
                       intervention_type=intervention_type)


        # FIX ME::: LESS ERROR IF EVERTHING IS AUTHOMATIC
        ## then we have to marginalize respect to the stochastic function (we have A=0 first and A=1 afterwards in the intervened data)
        W<-W0*(1-x$stochastic_probs[[protocol_name]][,match(treatment_variables[[j]],
                                                                   colnames(x$cumulative_stochastic_probs[[protocol_name]]))])+
               W1*x$stochastic_probs[[protocol_name]][,match(treatment_variables[[j]],
                                                                     colnames(x$cumulative_stochastic_probs[[protocol_name]]))]

      }

      else{
        # construction of clever covariates
        Wold <- rep(NA,length(Y))
        Wold[outcome_free_and_uncensored] <- lava::logit(fit_last)
      if (inherits(try(
        ## HERE
        # here if I understood correctly, in W we have the updated outcome obtained
        #by the logistic regression with the clever covariate and the fit_last as offset
        # and W has the dimension of Y (full data) with NAs if previously censored or had outcome already
        W <- update_Q(Y = Y,
                      logitQ = Wold,
                      cum.g = x$cumulative_intervention_probs[[protocol_name]][,censoring_variables[[j]]],
                      uncensored_undeterministic = outcome_free_and_uncensored,
                      intervention.match = x$intervention_match[[protocol_name]][,intervention_table[time == j-1]$variable],
                      intervention_type=intervention_type)
        ## W <- update_Q(Y = Y[!is.na(Y)],
        ## logitQ = lava::logit(fit_last[!is.na(Y)]),
        ## cum.g = x$cumulative_intervention_probs[[protocol_name]][,censoring_variables[[j]]][outcome_free_and_uncensored][!is.na(Y)],
        ## uncensored_undeterministic = outcome_free_and_uncensored,
        ## intervention.match = x$intervention_match[[protocol_name]][,intervention_table[time == j-1]$variable])
      ),"try-error"))
      stop(paste0("Fluctuation model used in the TMLE update step failed",
                  " in the attempt to run function update_Q at time point: ",j))
      # FIX ME: we have to take the weight at the correct time
      }
    }

    else{
      W <- fit_last
    }



    # FIXME: we do not need to save the following
    if (FALSE){
      x$sequential_outcome_regression[[target_name]]$predicted_values <- cbind(x$sequential_outcome_regression[[target_name]]$predicted_values,W)
      x$sequential_outcome_regression[[target_name]]$fit <- c(x$sequential_outcome_regression[[target_name]]$fit,list(fit_last))
      x$sequential_outcome_regression[[target_name]]$intervened_data <- c(x$sequential_outcome_regression[[target_name]]$intervened_data,list(intervened_data))
    }

    # calculate contribution to influence function
    # in case of the stochastic intervention we have in cumulative_intervention_probs already saved the weights
    if(intervention_type=="stochastic"){h.g.ratio<- x$cumulative_intervention_probs[[protocol_name]][,match(censoring_variables[[j]],colnames(x$cumulative_stochastic_probs[[protocol_name]]))] }
    else{
      h.g.ratio <- 1/x$cumulative_intervention_probs[[protocol_name]][,match(censoring_variables[[j]],colnames(x$cumulative_intervention_probs[[protocol_name]]))]
    }

    if (any(h.g.ratio>10000)) h.g.ratio <- pmin(h.g.ratio,10000)
    current_cnode <- as.character(x$prepared_data[[paste0(x$names$censoring,"_",j)]])

    index <- (current_cnode %in% x$names$uncensored_label) &  intervention_match[,intervention_table[time == j-1]$variable]
    if (any(h.g.ratio[index] != 0)) {
      x$IC[[target_name]][[protocol_name]][[label_time_horizon]][index] <- x$IC[[target_name]][[protocol_name]][[label_time_horizon]][index] + (Y[index] - W[index]) * h.g.ratio[index]
    }
    ## curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio,
    ## uncensored, intervention.match, regimes.with.positive.weight)
    # prepare next iteration
    set(x = x$prepared_data,
        i = which(outcome_free_and_uncensored), # atrisk: not censored, event-free, no competing risk
        j = "rtmle_predicted_outcome",
        value = W[which(outcome_free_and_uncensored)])
    # for those who have had an event or died or censored earlier
    set(x = x$prepared_data,
        i = which(!(outcome_free_and_uncensored)), # not atrisk
        j = "rtmle_predicted_outcome",
        # fixme j or j-1?
        value = x$prepared_data[[paste0(x$names$outcome,"_",j)]][which(!(outcome_free_and_uncensored))])
  }
  # g-formula and tmle estimator (since we do already before the multiplication with the stochastic function, when we do the mean it is not needed: CHRCK THIS)
  x$estimate[[target_name]][[protocol_name]][Time_horizon == time_horizon, Estimate := mean(x$prepared_data$rtmle_predicted_outcome)]
  ic <- x$IC[[target_name]][[protocol_name]][[label_time_horizon]] + x$prepared_data$rtmle_predicted_outcome - mean(x$prepared_data$rtmle_predicted_outcome)
  x$estimate[[target_name]][[protocol_name]][Time_horizon == time_horizon, Standard_error := sqrt(stats::var(ic)/NROW(x$prepared_data))]
  x$estimate[[target_name]][[protocol_name]][Time_horizon == time_horizon, Lower := Estimate-stats::qnorm(.975)*Standard_error]
  x$estimate[[target_name]][[protocol_name]][Time_horizon == time_horizon, Upper := Estimate+stats::qnorm(.975)*Standard_error]
  x$IC[[target_name]][[protocol_name]][[label_time_horizon]] <- ic
  return(x[])
}


######################################################################
### sequential_regression.R ends here
