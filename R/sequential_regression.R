### sequential_regression.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 30 2024 (14:30) 
## Version: 
## Last-Updated: Nov  1 2024 (07:35) 
##           By: Thomas Alexander Gerds
##     Update #: 139
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
sequential_regression <- function(x,
                                  target_name,
                                  protocol_name,
                                  time_horizon,
                                  learner,...){
    time = Time_horizon = Estimate = Standard_error = Lower = Upper = NULL
    # FIXME: inconsistent listing:
    intervention_table <- x$protocols[[protocol_name]]$intervention_table
    intervention_match <- x$intervention_match[[protocol_name]]
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
    for (j in reverse_time_scale){
        # subjects that are outcome-free and uncensored at the BEGINNING of the interval are 'atrisk'
        outcome_free_and_uncensored <- (x$followup$last_interval >= (j-1))
        # Note that at the time_horizon subjects who are censored in last interval have an NA outcome
        # but at any subsequent (lower time) we can predict the outcome also for subjects who are censored
        # during the current interval. So, while their observed outcome is NA their predicted_outcome is available
        # we keep the NA values and use outcome_free_and_uncensored rather than outcome_free_and_still_uncensored
        # because the TMLE update step needs to subset the predicted values to the non-NA outcomes
        if (j == time_horizon) {
            outcome_name <- outcome_variables[[time_horizon]]
            # FIXME: protocols could share the outcome formula?
            interval_outcome_formula = x$models[[protocol_name]][[outcome_variables[[j]]]]$formula
        } else {
            outcome_name <- "rtmle_predicted_outcome" 
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
        intervened_data = intervene(
            ## HERE
            ## data = x$prepared_data[outcome_free_and_uncensored][,1:(-1+match(outcome_variables[[j]],names(x$prepared_data))),with = FALSE],
            data = x$prepared_data[,intervenable_history,with = FALSE],
            intervention_table = intervention_table,
            time = j)
        # fit outcome regression
        # we always have the censoring variable in time interval _j *before* the outcome variable in same interval
        # hence to analyse Y_j we need C_j = "uncensored" for the modeling of Y_j
        # we thus remove subjects who are at risk at the beginning of the interval
        # but get censored during the interval (in fact, their outcome is NA)
        if (j == time_horizon)
            learn_variables <- c(history_of_variables,outcome_variables[[j]])
        else
            learn_variables <- c(history_of_variables,"rtmle_predicted_outcome")
        
        args <- list(character_formula = interval_outcome_formula,
                     data = x$prepared_data[outcome_free_and_uncensored][,learn_variables,with = FALSE],
                     intervened_data = intervened_data[outcome_free_and_uncensored],...)
        if (length(learner)>1){
            if (j == time_horizon)
                args <- c(args,list(learners = learner, outcome_variable = outcome_variables[[j]], id_variable = x$names$id))
            else
                args <- c(args,list(learners = learner, outcome_variable = "rtmle_predicted_outcome", id_variable = x$names$id))
            if (inherits(try(
                fit_last <- do.call("superlearn",args),silent = FALSE),
                "try-error")) {
                ## browser(skipCalls=1L)
                stop(paste0("Sequential regression fit failed with formula:\n",interval_outcome_formula))
            }
        }else{
            if (inherits(try(
                fit_last <- do.call(learner,args),silent = FALSE),
                "try-error")) {
                ## browser(skipCalls=1L)
                stop(paste0("Sequential regression fit failed with formula:\n",interval_outcome_formula))
            }
        }
        # save fitted object
        x$models[[protocol_name]][[outcome_variables[[j]]]]$fit <- attr(fit_last,"fit")
        data.table::setattr(fit_last,"fit",NULL)
        fit_last <- as.numeric(fit_last)
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
            ## current_prediction <- lava::logit(fit_last[!is.na(Y)])
            ## current_outcome <- Y[!is.na(Y)]
            ## print(data.table(Y = Y,fit_last = fit_last))
            # construction of clever covariates
            Wold <- rep(NA,length(Y))
            Wold[outcome_free_and_uncensored] <- lava::logit(fit_last)
            ## if (j == 1) browser(skipCalls=1L)
            if (inherits(try(
                ## HERE
                W <- update_Q(Y = Y,
                              logitQ = Wold,
                              cum.g = x$cumulative_intervention_probs[[protocol_name]][,censoring_variables[[j]]],
                              uncensored_undeterministic = outcome_free_and_uncensored,
                              intervention.match = x$intervention_match[[protocol_name]][,intervention_table[time == j-1]$variable])
                ## W <- update_Q(Y = Y[!is.na(Y)],
                ## logitQ = lava::logit(fit_last[!is.na(Y)]),
                ## cum.g = x$cumulative_intervention_probs[[protocol_name]][,censoring_variables[[j]]][outcome_free_and_uncensored][!is.na(Y)], 
                ## uncensored_undeterministic = outcome_free_and_uncensored,
                ## intervention.match = x$intervention_match[[protocol_name]][,intervention_table[time == j-1]$variable])
            ),"try-error"))
                stop(paste0("Fluctuation model used in the TMLE update step failed",
                            " in the attempt to run function update_Q at time point: ",j))
        }else{
            W <- fit_last
        }
        # FIXME: we do not need to save the following
        if (FALSE){
            x$sequential_outcome_regression[[target_name]]$predicted_values <- cbind(x$sequential_outcome_regression[[target_name]]$predicted_values,W)
            x$sequential_outcome_regression[[target_name]]$fit <- c(x$sequential_outcome_regression[[target_name]]$fit,list(fit_last))
            x$sequential_outcome_regression[[target_name]]$intervened_data <- c(x$sequential_outcome_regression[[target_name]]$intervened_data,list(intervened_data))
        }
        # calculate contribution to influence function
        h.g.ratio <- 1/x$cumulative_intervention_probs[[protocol_name]][,match(censoring_variables[[j]],colnames(x$cumulative_intervention_probs[[protocol_name]]))]
        if (any(h.g.ratio>10000)) h.g.ratio <- pmin(h.g.ratio,10000)
        current_cnode <- as.character(x$prepared_data[[paste0(x$names$censoring,"_",j)]])
        index <- (current_cnode%in%x$names$uncensored_label) & intervention_match[,intervention_table[time == j-1]$variable]
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
    # g-formula and tmle estimator
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
