### sequential_regression.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 30 2024 (14:30)
## Version:
## Last-Updated: feb 26 2026 (12:59) 
##           By: Thomas Alexander Gerds
##     Update #: 547
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
                                  learner,
                                  estimator,
                                  seed = seed){
    time = Target = Protocol = Time_horizon = Estimator = Estimate = Target_parameter = Standard_error = Lower = Upper = rtmle_predicted_outcome = NULL
    N <- NROW(x$prepared_data)
    intervention_table <- x$protocols[[protocol_name]]$intervention_table
    intervention_match <- x$protocols[[protocol_name]]$intervention_match
    if (length(x$names$censoring)>0){
        censoring_variables <- paste0(x$names$censoring,"_",1:time_horizon)
    }else{
        censoring_variables <- NULL
    }
    competing_variables <- paste0(x$names$competing,"_",1:(time_horizon-1))
    outcome_variables <- paste0(x$names$outcome,"_",1:time_horizon)
    # the first step is the outcome regression of the last time interval
    x$sequential_outcome_regression[[target_name]] = vector(3,mode = "list")
    stopifnot(time_horizon>0)
    label_time_horizon <- paste0("time_horizon_",time_horizon)
    reverse_time_scale <- rev(seq(1,time_horizon,1))
    for (k in reverse_time_scale){
        # subjects that are outcome-free and uncensored at the BEGINNING of the interval are 'atrisk'
        if (length(x$followup) == 0){
            outcome_free_and_uncensored <- rep(TRUE,N)
        }else{
            outcome_free_and_uncensored <- (x$followup$last_interval >= (k-1))
        }
        # we always have the censoring variable in time interval _k *before* the outcome variable in same interval
        # hence to analyse Y_k we need C_k = "uncensored" for the modeling of Y_k
        # we thus remove subjects who are at risk at the beginning of the interval
        # but get censored during the interval (in fact, their outcome is NA)
        # without covariates what happens is that the mean outcome probability is imputed
        # for those who are censored during this interval.
        if (length(x$names$censoring)>0){
            current_cnode <- as.character(x$prepared_data[[paste0(x$names$censoring,"_",k)]])
            outcome_free_and_uncensored_outcome <- outcome_free_and_uncensored & (current_cnode %in% x$names$uncensored_label)
        }else{
            outcome_free_and_uncensored_outcome <- outcome_free_and_uncensored
        }
        # Note that at the time_horizon subjects who are censored in last interval have an NA outcome
        # but at any subsequent step of the algorithm, that is in an earlier time interval,
        # the predicted outcome is available also for subjects who are censored
        # during the current interval. So, while their observed outcome is NA
        # their predicted_outcome is available.
        interval_outcome_formula <- x$models[[paste0("time_",(k-1))]][["outcome"]][[1]]$formula
        if (k == time_horizon) {        
            outcome_name <- outcome_variables[[time_horizon]]
        }else{
            outcome_name <- "rtmle_predicted_outcome"
            interval_outcome_formula <- stats::update(stats::formula(interval_outcome_formula),"rtmle_predicted_outcome~.")
        }
        Y <- x$prepared_data[[outcome_name]]
        # intervene according to protocol for targets
        intervened_data <- do.call(x$protocol[[protocol_name]]$intervene_function,
                                   list(data = x$prepared_data[outcome_free_and_uncensored],
                                        intervention_table = intervention_table,
                                        time = k))
        # fit outcome regression
        fit_last_interval <- fitter(intervention_node = k,
                                    learner = learner,
                                    formula = interval_outcome_formula,
                                    data = x$prepared_data[outcome_free_and_uncensored_outcome],
                                    intervened_data = intervened_data,
                                    id_variable = x$names$id,
                                    minority_threshold = x$tuning_parameters$minority_threshold,
                                    seed = seed,
                                    diagnostics = x$diagnostics)
        # save fitted object
        x$models[[paste0("time_",(k-1))]][["outcome"]][[protocol_name]]$fit <- attr(fit_last_interval,"fit")
        x$diagnostics <- attr(fit_last_interval,"diagnostics")
        # remove attributes
        fit_last_interval <- as.numeric(fit_last_interval)
        # set predicted value as outcome for next regression
        if (estimator == "tmle"){
            old_action <- options()$na.action
            on.exit(options(na.action = old_action))
            options(na.action = "na.pass")
            # note that we can still predict those who are censored at C_j but uncensored at C_{j-1}
            ## y <- riskRegression::predictRisk(fit_last_interval,newdata = intervened_data)
            # avoid missing values due to logit
            if (any(fit_last_interval[!is.na(fit_last_interval)] <= 0)) fit_last_interval <- pmax(fit_last_interval,x$tuning_parameters$prediction_range[1])
            if (any(fit_last_interval[!is.na(fit_last_interval)] >= 1)) fit_last_interval <- pmin(fit_last_interval,x$tuning_parameters$prediction_range[2])
            # TMLE update step
            ipos <- x$protocols[[protocol_name]]$intervention_last_nodes[k]
            inverse_probability_weights <- x$protocols[[protocol_name]]$cumulative_intervention_probs[,ipos]
            # weight truncation
            if (is.numeric(x$tuning_parameters$weight_truncation)){
                inverse_probability_weights <- pmax(pmin(inverse_probability_weights,
                                                         x$tuning_parameters$weight_truncation[2]),
                                                    x$tuning_parameters$weight_truncation[1])
            }
            # the column names A_1,B_1,E_1 of the intervention_match table are made with paste
            # in function intervention_probabilities
            intervention_node_name <- paste(intervention_table[time == k-1]$variable,collapse = ",")
            # use only data from subjects who are uncensored in current interval
            # construction of clever covariates
            predicted_outcome_previous <- rep(NA,length(Y))
            predicted_outcome_previous[outcome_free_and_uncensored] <- stats::qlogis(fit_last_interval)
            if (nchar(intervention_node_name)>0){
                imatch <- (intervention_match[,intervention_node_name]%in% 1)
            }else{
                imatch <- rep(1,N)
                imatch[!outcome_free_and_uncensored] <- NA
            }
            if (any(inverse_probability_weights[!is.na(Y) & outcome_free_and_uncensored & as.vector(imatch)] == 0)){
                stop("Exactly zero intervention probabilities encountered at the attempt to run the TMLE-update fluctuation model.\nYou may have to consider changing the target parameter or bounding the intervention probabilities somehow.\nGood luck!")
            }
            if (inherits(try(
                predicted_outcome <- tmle_update(Y = Y,
                                                 offset = predicted_outcome_previous,
                                                 intervention_probs = inverse_probability_weights,
                                                 outcome_free_and_uncensored = outcome_free_and_uncensored,
                                                 intervention_match = imatch,
                                                 k = k,protocol = protocol_name)
            ),"try-error")){
                stop(paste0("Fluctuation model used in the TMLE update step faile",
                            " in the attempt to run function tmle_update at time point: ",k))
            }
            ## FIXME: why do we NOT need the following?
            ## predicted_outcome[!outcome_free_and_uncensored] <- Y[!outcome_free_and_uncensored]
            #
            # calculate contribution to influence function
            #
            IPW <- 1/inverse_probability_weights
            # FIXME: what is this constant 10000? 
            ## if (any(IPW>10000)) IPW <- pmin(IPW,10000)
            ## if (max(IPW,na.rm = TRUE)>nrow(x$prepared_data)) IPW <- pmin(IPW,nrow(x$prepared_data))
            if (length(x$names$censoring)>0){
                if (nchar(intervention_node_name)>0){
                    index <- (current_cnode%in%x$names$uncensored_label) & (intervention_match[,intervention_node_name] %in% 1)
                }else{
                    index <- (current_cnode%in%x$names$uncensored_label)
                }
            }else{
                if (nchar(intervention_node_name)>0){
                    index <- (intervention_match[,intervention_node_name] %in% 1)
                }else{
                    index <- rep(1,N)
                    index[!outcome_free_and_uncensored] <- NA
                }
            }
            if (any(IPW[index] != 0)) {
                x$IC[[target_name]][[protocol_name]][[label_time_horizon]][index] <- x$IC[[target_name]][[protocol_name]][[label_time_horizon]][index] + (Y[index] - predicted_outcome[index]) * IPW[index]
            }
        } else {
            # g-formula
            predicted_outcome <- x$prepared_data[[paste0(x$names$outcome,"_",k)]]
            predicted_outcome[which(outcome_free_and_uncensored)] <- fit_last_interval
            # FIXME: how to estimate the influence function without inverse probability weights?
        }
        # prepare next iteration
        set(x = x$prepared_data,
            i = which(outcome_free_and_uncensored), # atrisk: not censored, event-free, no competing risk
            j = "rtmle_predicted_outcome",
            value = predicted_outcome[which(outcome_free_and_uncensored)])
        # for those who have had an event or died or censored earlier
        # the previous value is set
        set(x = x$prepared_data,
            i = which(!(outcome_free_and_uncensored)), # not atrisk
            j = "rtmle_predicted_outcome",
            value = x$prepared_data[[paste0(x$names$outcome,"_",k)]][which(!(outcome_free_and_uncensored))])
    }
    target_parameter <- "Risk"
    # g-formula and tmle estimator
    ic <- x$IC[[target_name]][[protocol_name]][[label_time_horizon]] + x$prepared_data$rtmle_predicted_outcome - mean(x$prepared_data$rtmle_predicted_outcome)
    SE = sqrt(stats::var(ic)/N)
    x$estimate[["Main_analysis"]][Target == target_name & Protocol ==  protocol_name & Time_horizon == time_horizon & Target_parameter == target_parameter & Estimator == estimator,
                                  `:=`(Estimate = mean(x$prepared_data$rtmle_predicted_outcome),
                                       Standard_error = SE)]
    x$estimate[["Main_analysis"]][Target == target_name & Protocol ==  protocol_name & Time_horizon == time_horizon & Target_parameter == target_parameter & Estimator == estimator,
                                  `:=`(Lower = Estimate-stats::qnorm(.975)*SE,
                                       Upper = Estimate+stats::qnorm(.975)*SE)][]
    x$IC[[target_name]][[protocol_name]][[label_time_horizon]] <- ic
    # clean up for the next run
    data.table::setkey(x$estimate$Main_analysis,Target,Protocol,Target_parameter,Time_horizon,Estimator)
    x$prepared_data[,rtmle_predicted_outcome := NULL]
    # NOTE if we would return x[] instead of x then x looses its class!
    return(x)
}


######################################################################
### sequential_regression.R ends here
