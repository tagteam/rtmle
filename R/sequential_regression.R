### sequential_regression.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 30 2024 (14:30) 
## Version: 
## Last-Updated: Oct 16 2024 (15:58) 
##           By: Thomas Alexander Gerds
##     Update #: 46
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
        # formula
        # FIXME: protocols could share the outcome formula?
        interval_outcome_formula = stats::formula(x$models[[protocol_name]][["outcome"]][[outcome_variables[[j]]]]$formula)
        #  at the last time interval the observed outcome is used
        if (j != time_horizon) {
            interval_outcome_formula <- stats::update(interval_outcome_formula,"predicted_outcome~.")
            Yhat <- x$prepared_data[["predicted_outcome"]]
        }else{
            Yhat <- x$prepared_data[[all.vars(interval_outcome_formula)[[1]]]]
        }
        # subjects that are outcome-free and uncensored at the beginning of the interval are 'atrisk'
        # subjects not censored within interval are 'uncensored_interval'
        # note that subjects who are uncensored_interval are also atrisk
        # because the prepare_data function has set the values of the censoring
        # variables that are later than an outcome or competing risk event to NA
        if (j > 1){
            outcome_free <- x$prepared_data[[outcome_variables[[j-1]]]]%in%0
            # competing risks
            if (length(competing_variables)>0){
                outcome_free <- outcome_free&x$prepared_data[[competing_variables[[j-1]]]]%in%0
            }
            uncensored <- x$prepared_data[[censoring_variables[[j-1]]]]%in%x$names$uncensored_label
        }else{
            outcome_free <- rep(TRUE,nrow(x$prepared_data))
            uncensored <- rep(TRUE,nrow(x$prepared_data))
        }
        # fit outcome regression
        # we always have the censoring variable in time interval _j *before* the outcome variable in same interval
        # hence to analyse Y_j we need C_j = x$names$uncensored_label for the modeling of Y_j
        # we thus remove subjects who are at risk at the beginning of the interval
        # but then get censored because their outcome is NA
        fit_last <- do.call(learner,list(formula = interval_outcome_formula,data = x$prepared_data[outcome_free&x$prepared_data[[censoring_variables[[j]]]]%in%x$names$uncensored_label],...))
        # save fitted object
        x$models[[protocol_name]][["outcome"]][[outcome_variables[[j]]]]$fit <- fit_last
        # intervene according to protocol for targets
        # FIXME: intervene on all variables or only those after
        #        time j? those in current outcome_formula
        intervened_data = intervene(
            ## formula = interval_outcome_formula,
            ## data = x$prepared_data[outcome_free&uncensored],
            data = x$prepared_data,
            intervention_table = intervention_table,
            time = j)
        # set predicted value as outcome for next regression
        ## FIXME: if (target$estimator == "tmle")
        old_action <- options()$na.action
        on.exit(options(na.action = old_action))
        options(na.action = "na.pass")
        # note that we can still predict those who are censored at C_j but uncensored at C_{j-1}
        y <- riskRegression::predictRisk(fit_last,newdata = intervened_data)
        # avoid missing values due to logit
        if (x$targets[[target_name]]$estimator == "tmle"){
            if (any(y[!is.na(y)] < 0)) y <- pmax(y,0.0001)
            if (any(y[!is.na(y)] > 1)) y <- pmin(y,0.9999)
        }
        ## y <- pmax(pmin(y,0.99999),0.00001)
        ## y <- pmax(pmin(y,0.9999),0.0001)
        # FIXME: decide to either use paste0(x$names$censoring,"_",j) or censoring_variables[[j]]
        current_cnode = as.character(x$prepared_data[[paste0(x$names$censoring,"_",j)]])
        if (length(x$targets[[target_name]]$estimator) == 0 || x$targets[[target_name]]$estimator == "tmle"){
            if (inherits(try(W <- update_Q(Y = Yhat,
                                           logitQ = lava::logit(y),
                                           cum.g = x$cumulative_intervention_probs[[protocol_name]][,censoring_variables[[j]]], 
                                           uncensored_undeterministic = outcome_free & (current_cnode%in%x$names$uncensored_label),
                                           intervention.match = x$intervention_match[[protocol_name]][,intervention_table[time == j-1]$variable])),"try-error"))
                stop("Fluctuation model used in the TMLE update step failed in the attempt to run function update_Q")
        }else{
            W <- y
        }
        x$sequential_outcome_regression[[target_name]]$predicted_values <- cbind(x$sequential_outcome_regression[[target_name]]$predicted_values,W)
        x$sequential_outcome_regression[[target_name]]$fit <- c(x$sequential_outcome_regression[[target_name]]$fit,list(fit_last))
        x$sequential_outcome_regression[[target_name]]$intervened_data <- c(x$sequential_outcome_regression[[target_name]]$intervened_data,list(intervened_data))
        # calculate contribution to influence function
        h.g.ratio <- 1/x$cumulative_intervention_probs[[protocol_name]][,match(censoring_variables[[j]],colnames(x$cumulative_intervention_probs[[protocol_name]]))]
        index <- (current_cnode%in%x$names$uncensored_label) & intervention_match[,intervention_table[time == j-1]$variable]
        if (any(h.g.ratio[index] != 0)) {
            x$IC[[target_name]][[protocol_name]][[label_time_horizon]][index] <- x$IC[[target_name]][[protocol_name]][[label_time_horizon]][index] + (Yhat[index] - W[index]) * h.g.ratio[index]
        }
        ## curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio,
        ## uncensored, intervention.match, regimes.with.positive.weight)
        # prepare next iteration
        set(x = x$prepared_data,
            i = which(outcome_free&uncensored), # atrisk: not censored, event-free, no competing risk
            j = "predicted_outcome",
            value = W[which(outcome_free&uncensored)])
        # for those who have had an event or died or censored earlier
        set(x = x$prepared_data,
            i = which(!(outcome_free&uncensored)), # not atrisk
            j = "predicted_outcome",
            # fixme j or j-1?
            value = x$prepared_data[[paste0(x$names$outcome,"_",j)]][which(!(outcome_free&uncensored))])
    }
    ## if (time_horizon == 2) browser()
    # g-formula and tmle estimator
    x$estimate[[target_name]][[protocol_name]][Time_horizon == time_horizon, Estimate := mean(x$prepared_data$predicted_outcome)]
    ic <- x$IC[[target_name]][[protocol_name]][[label_time_horizon]] + x$prepared_data$predicted_outcome - mean(x$prepared_data$predicted_outcome)
    x$estimate[[target_name]][[protocol_name]][Time_horizon == time_horizon, Standard_error := sqrt(stats::var(ic)/NROW(x$prepared_data))]
    x$estimate[[target_name]][[protocol_name]][Time_horizon == time_horizon, Lower := Estimate-stats::qnorm(.975)*Standard_error]
    x$estimate[[target_name]][[protocol_name]][Time_horizon == time_horizon, Upper := Estimate+stats::qnorm(.975)*Standard_error]
    x$IC[[target_name]][[protocol_name]][[label_time_horizon]] <- ic
    return(x[])
}


######################################################################
### sequential_regression.R ends here
