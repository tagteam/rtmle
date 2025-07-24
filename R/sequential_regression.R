### sequential_regression.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 30 2024 (14:30)
## Version:
## Last-Updated: Jul 24 2025 (11:41) 
##           By: Thomas Alexander Gerds
##     Update #: 338
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
                                  seed = seed,
                                  ...){
    time = Target = Protocol = Time_horizon = Estimate = Target_parameter = Standard_error = Lower = Upper = NULL
    N <- NROW(x$prepared_data)
    # FIXME: inconsistent listing:
    intervention_table <- x$protocols[[protocol_name]]$intervention_table
    intervention_match <- x$intervention_match[[protocol_name]]
    # FIXME: treatment_variables=intervention_table$variable?
    treatment_variables <- x$protocols[[protocol_name]]$intervention_table[time%in%0:time_horizon]$variable
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
    for (j in reverse_time_scale){
        # subjects that are outcome-free and uncensored at the BEGINNING of the interval are 'atrisk'
        if (j == 1){
            outcome_free_and_uncensored <- rep(TRUE,N)
        }else{
            outcome_free_and_uncensored <- (x$followup$last_interval >= (j-1))
        }
        # Note that at the time_horizon subjects who are censored in last interval have an NA outcome
        # but at any subsequent (lower time) we can predict the outcome also for subjects who are censored
        # during the current interval. So, while their observed outcome is NA their predicted_outcome is available
        # we keep the NA values and use outcome_free_and_uncensored rather than outcome_free_and_still_uncensored
        # because the TMLE update step needs to subset the predicted values to the non-NA outcomes
        outcome_name <- outcome_variables[[time_horizon]]
        # FIXME: delete the if query when obsolete way of specifying formulas is gone
        if (protocol_name %in% names(x$models)){
            interval_outcome_formula = x$models[[protocol_name]][[outcome_variables[[j]]]]$formula
        }else{
            interval_outcome_formula = x$models[[outcome_variables[[j]]]]$formula
        }
        if (j < time_horizon) {
            outcome_name <- "rtmle_predicted_outcome"
            interval_outcome_formula <- stats::update(stats::formula(interval_outcome_formula),"rtmle_predicted_outcome~.")
        }
        ## HERE
        ## Y <- x$prepared_data[outcome_free_and_uncensored][[outcome_name]]
        Y <- x$prepared_data[[outcome_name]]
        # intervene according to protocol for targets
        # FIXME: intervene on all variables or only those after
        #        time j? those in current outcome_formula
        # FIXME: remove id variable below here
        history_of_variables <- names(x$prepared_data)[1:(-1+match(outcome_variables[[j]],names(x$prepared_data)))]
        intervenable_history <- setdiff(history_of_variables,c(outcome_variables,censoring_variables,competing_variables))
        intervened_data <- do.call(x$protocol[[protocol_name]]$intervene_function,
                                   list(data = x$prepared_data[,intervenable_history,with = FALSE],
                                        intervention_table = intervention_table,
                                        time = j))
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
        if (NROW(current_data) == 0) {
            stop("No data available for g-estimation")
        }
        
        current_constants <- sapply(current_data, function(x){length(unique(x))==1})
        if (any(current_constants)) {
            current_constants <- names(current_constants[current_constants])
        }else{
            current_constants <- NULL
        }
        # check if there is variability in the outcome variable
        if (outcome_variables[[j]]%in%current_constants)
            stop(paste0("No variation in the outcome variable: ",
                        outcome_variables[[j]],
                        ".\nCheck the prepared data:\nx$prepared_data[x$followup$last_interval>=",
                        j-1,",.(",x$names$id,",",outcome_variables[[j]],")]"))
        # remove constant predictor variables
        interval_outcome_formula_vars <- all.vars(stats::formula(interval_outcome_formula))
        if (length(current_constants)>0){
            # FIXME: table warnings and show number of times per warning
            ## x$warnings <- paste0("Removing constant variables at time ",j,":\n",paste0(current_constants,collapse = ", "))
            interval_outcome_formula <- delete_variables_from_formula(character_formula = interval_outcome_formula,delete_vars = current_constants)
        }
        args <- list(character_formula = interval_outcome_formula,
                     data = current_data[,!(names(current_data)%in%current_constants),with = FALSE],
                     intervened_data = intervened_data[outcome_free_and_uncensored],...)
        # super learner needs name of outcome variable
        if (length(learner)>1){
            if (j == time_horizon){
                args <- c(args,list(learners = learner,
                                    outcome_variable = outcome_variables[[j]],
                                    outcome_target_level = NULL,
                                    id_variable = x$names$id))
            }else{
                args <- c(args,list(learners = learner,
                                    outcome_variable = "rtmle_predicted_outcome",
                                    outcome_target_level = NULL,
                                    id_variable = x$names$id))
            }
            if (inherits(try(
                fit_last <- do.call("superlearn",c(args,list(seed = seed))),silent = FALSE),
                "try-error")) {
                stop(paste0("Sequential regression fit failed with formula:\n",interval_outcome_formula))
            }
        }else{
            # take care of case where additional arguments are passed to a single learner
            if (is.list(learner)){
                learner_args <- learner[[1]][names(learner[[1]]) != "learner_fun"]
                args <- c(args,learner_args)
                learner_fun <- learner[[1]][["learner_fun"]]
            }else{
                learner_fun <- learner
            }
            # single learners do not need the name of outcome variable
            if (inherits(try(
                fit_last <- do.call(learner_fun,args),silent = FALSE),
                "try-error")) {
                stop(paste0("Sequential regression fit failed with formula:\n",interval_outcome_formula))
            }
        }
        # save fitted object
        x$models[[outcome_variables[[j]]]]$fit <- attr(fit_last,"fit")
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
            if (any(fit_last[!is.na(fit_last)] >= 1)) fit_last <- pmin(fit_last,0.9999)
        }
        ## y <- pmax(pmin(y,0.99999),0.00001)
        ## y <- pmax(pmin(y,0.9999),0.0001)
        # TMLE update step
        if (length(x$names$censoring)>0){
            ipos <- censoring_variables[[j]]
        }else{
            # FIXME: this is not easy to read, but when there are multiple treatment variables
            #        at time j, such as A_j and B_j then we match on the last
            ipos <- rev(treatment_variables[[j]])[[1]]
        }
        # the column names A_1,B_1,E_1 of the intervention_match table are made with paste
        # in function intervention_probabilities
        intervention_node_name <- paste(intervention_table[time == j-1]$variable,collapse = ",")
        if (length(x$targets[[target_name]]$estimator) == 0 || x$targets[[target_name]]$estimator == "tmle"){
            # use only data from subjects who are uncensored in current interval
            # construction of clever covariates
            W_previous <- rep(NA,length(Y))
            W_previous[outcome_free_and_uncensored] <- lava::logit(fit_last)
            # FIXME: this test of infinite inverse_probability_weights needs more work
            inverse_probability_weights <- x$cumulative_intervention_probs[[protocol_name]][,ipos]
            imatch <- (x$intervention_match[[protocol_name]][,intervention_node_name]%in% 1)
            if (any(inverse_probability_weights[!is.na(Y) & outcome_free_and_uncensored & as.vector(imatch)] == 0)){
                stop("Exactly zero intervention probabilities encountered at the attempt to run the TMLE-update fluctuation model.\nYou may have to consider changing the target parameter or bounding the intervention probabilities somehow.\nGood luck!")
            }
            if (inherits(try(
                W <- tmle_update(Y = Y,
                                 offset = W_previous,
                                 intervention_probs = inverse_probability_weights,
                                 outcome_free_and_uncensored = outcome_free_and_uncensored,
                                 intervention_match = imatch)
            ),"try-error"))
                stop(paste0("Fluctuation model used in the TMLE update step failed",
                            " in the attempt to run function tmle_update at time point: ",j))
            ## FIXME: why do we not need the following?
            ## W[!outcome_free_and_uncensored] <- Y[!outcome_free_and_uncensored]
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
        h.g.ratio <- 1/x$cumulative_intervention_probs[[protocol_name]][,ipos]
        if (any(h.g.ratio>10000)) h.g.ratio <- pmin(h.g.ratio,10000)
        if (length(x$names$censoring)>0){
            current_cnode <- as.character(x$prepared_data[[paste0(x$names$censoring,"_",j)]])
            index <- (current_cnode%in%x$names$uncensored_label) & (intervention_match[,intervention_node_name] %in% 1)
        }else{
            index <- (intervention_match[,intervention_node_name] %in% 1)
        }
        
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
    target_parameter <- "Risk"
    # g-formula and tmle estimator
    ic <- x$IC[[target_name]][[protocol_name]][[label_time_horizon]] + x$prepared_data$rtmle_predicted_outcome - mean(x$prepared_data$rtmle_predicted_outcome)
    SE = sqrt(stats::var(ic)/N)
    x$estimate[["Main_analysis"]][Target == target_name & Protocol ==  protocol_name & Time_horizon == time_horizon & Target_parameter == target_parameter,
                                  Estimate := mean(x$prepared_data$rtmle_predicted_outcome)]
    x$estimate[["Main_analysis"]][Target == target_name & Protocol ==  protocol_name & Time_horizon == time_horizon & Target_parameter == target_parameter,
                                  `:=`(Standard_error = SE,
                                       Lower = Estimate-stats::qnorm(.975)*SE,
                                       Upper = Estimate+stats::qnorm(.975)*SE)]
    # FIXME: is there a better way to circumvent data.tables
    #        print after := problem?
    x$estimate[["Main_analysis"]][]
    x$IC[[target_name]][[protocol_name]][[label_time_horizon]] <- ic
    # NOTE if we would return x[] instead of x then x looses its class!
    return(x)
}


######################################################################
### sequential_regression.R ends here
