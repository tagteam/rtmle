### run_rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  1 2024 (09:11) 
## Version: 
## Last-Updated: Aug  1 2024 (11:13) 
##           By: Thomas Alexander Gerds
##     Update #: 201
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @export
run_rtmle <- function(x,targets,refit = FALSE,speed = FALSE,...){
    # FIXME: for target in targets
    available_targets <- names(x$targets)
    if (!missing(targets)) {
        if (!(all(targets %in% available_targets)))
            stop(paste0("Requested targets: \n",
                        paste(targets,collapse = ","),
                        "\n\navailable targets:\n",
                        paste(available_targets,collapse = ",")))
        run_these_targets <- intersect(targets,available_targets)
    }else{
        run_these_targets <- available_targets
    }
    for (target_name in run_these_targets){
        message("Running target: ",target_name)
        ## protocols <- names(x$protocols)
        protocols <- x$targets[[target_name]]$protocols
        for (protocol_name in protocols){
            message("Current protocol: ",protocol_name)
            protocol_data <- x$prepared_data[[protocol_name]]$data
            # initialize influence curve vector
            x$IC[[target_name]][[protocol_name]] <- numeric(nrow(protocol_data))
            #
            # g-part: fit nuisance parameter models for propensity and censoring
            #
            current_protocol <- x$protocols[[protocol_name]]
            intervention_table <- current_protocol$intervention_table
            protocol_Anodes <- x$prepared_data[[protocol_name]]$Anodes
            protocol_Cnodes <- x$prepared_data[[protocol_name]]$Cnodes
            protocol_Dnodes <- x$prepared_data[[protocol_name]]$Dnodes
            protocol_Ynodes <- x$prepared_data[[protocol_name]]$Ynodes
            #
            # define a matrix with the cumulative intervention/censoring probabilities
            #
            if (refit || length(x$cumulative_intervention_probs[[protocol_name]]) == 0){
                intervention_probs <- data.table(ID = protocol_data[[x$name_id]])
                intervention_probs <- cbind(intervention_probs,matrix(1,nrow = NROW(protocol_data),ncol = length(protocol_Anodes)+length(protocol_Cnodes)))
                setnames(intervention_probs,c(x$name_id,c(rbind(protocol_Anodes,protocol_Cnodes))))
                # predict the propensity score/1-probability of censored
                # intervene according to protocol for targets
                # in the last time interval we do not need propensities/censoring probabilities  
                for (j in x$times[-length(x$times)]){
                    if (j > 0){
                        outcome_free <- protocol_data[[protocol_Ynodes[[j]]]]%in%0
                        # competing risks
                        if (length(protocol_Dnodes)>0) outcome_free <- outcome_free&protocol_data[[protocol_Dnodes[[j]]]]%in%0
                        uncensored <- protocol_data[[protocol_Cnodes[[j]]]]%in%"uncensored"
                    }else{
                        # assume that all are at risk at time 0
                        outcome_free <- rep(TRUE,nrow(protocol_data))
                        uncensored <- rep(TRUE,nrow(protocol_data))
                    }
                    if (any(outcome_free&uncensored)){
                        # fit the propensity regression model
                        if (refit || length(x$models[[protocol_name]][["propensity"]][[protocol_Anodes[[j+1]]]]$fit) == 0){
                            if (!is.null(ff <- x$models[[protocol_name]][["propensity"]][[protocol_Anodes[[j+1]]]]$formula)){
                                x$models[[protocol_name]][["propensity"]][[protocol_Anodes[[j+1]]]]$fit <- fit_propensity_models(formula = ff,
                                                                                                                                 data = protocol_data[outcome_free&uncensored],
                                                                                                                                 speed = speed)
                            }
                        }
                        # fit censoring model
                        if (refit || length(x$models[[protocol_name]][["censoring"]]$fit) == 0){
                            if (!is.null(ff <- x$models[[protocol_name]][["censoring"]][[protocol_Cnodes[[j+1]]]]$formula)){
                                x$models[[protocol_name]][["censoring"]][[protocol_Cnodes[[j+1]]]]$fit <- fit_censoring_models(formula = ff,
                                                                                                                               data = protocol_data[outcome_free&uncensored],
                                                                                                                               speed = speed)
                            }
                        }
                    }
                    # set the treatment variables to their protocol values
                    ## FIXME: could subset data to the variables in the current formula
                    # FIXME: also those censored?
                    intervened_data = intervene(data = protocol_data[outcome_free],
                                                intervention_table = intervention_table,
                                                time = j)
                    intervention_probs[outcome_free][[protocol_Anodes[[j+1]]]] <- predict(x$models[[protocol_name]][["propensity"]][[protocol_Anodes[[j+1]]]]$fit,
                                                                                          newdata = intervened_data,
                                                                                          type = "response")
                    intervention_probs[outcome_free][[protocol_Cnodes[[j+1]]]] <- predict(x$models[[protocol_name]][["censoring"]][[protocol_Cnodes[[j+1]]]]$fit,newdata = intervened_data,type = "response")
                }
                x$cumulative_intervention_probs[[protocol_name]] <- matrixStats::rowCumprods(as.matrix(intervention_probs[,-1,with = FALSE]))
            }
            #
            # define a matrix which indicates if the intervention is followed
            #
            if (length(x$intervention_match[[protocol_name]]) == 0){
                # fixme for target in targets
                intervention_match <- matrix(0,ncol = length(protocol_Anodes),nrow = nrow(protocol_data))
                for(j in x$times[-c(length(x$times))]){
                    if (j == 0)
                        intervention_match[,j+1] <- previous <- (protocol_data[[intervention_table[j+1]$variable]] %in% c(intervention_table[j+1]$value,NA))
                    else
                        intervention_match[,j+1] <- previous <- previous*(protocol_data[[intervention_table[j+1]$variable]] %in% c(intervention_table[j+1]$value,NA))
                }
                colnames(intervention_match) <- protocol_Anodes
                x$intervention_match[[protocol_name]] <- intervention_match
            }
            # 
            # Q-part: loop backwards in time through iterative condtional expectations
            #
            # the first step is the outcome regression of the last time interval
            for (j in rev(x$times)[-length(x$times)]){
                # formula
                # FIXME: protocols could share the outcome formula?
                interval_outcome_formula = formula(x$models[[protocol_name]][["outcome"]][[protocol_Ynodes[[j]]]]$formula)
                #  at the last time interval the observed outcome is used
                if (j != max(x$times)) {
                    interval_outcome_formula <- update(interval_outcome_formula,"predicted_outcome~.")
                    Yhat <- protocol_data[["predicted_outcome"]]
                }else{
                    Yhat <- protocol_data[[all.vars(interval_outcome_formula)[[1]]]]
                }
                # data at-risk at the beginning of the interval
                if (j > 1){
                    outcome_free <- protocol_data[[paste0(x$name_outcome,"_",j-1)]]%in%0
                    # competing risks
                    # FIXME: should be protocol_Dnodes[[j-1]] or [[j]]? instead of paste0
                    if (length(protocol_Dnodes)>0) outcome_free <- outcome_free&protocol_data[[paste0(x$name_competing,"_",j-1)]]%in%0
                    uncensored <- protocol_data[[paste0(x$name_censoring,"_",j-1)]]%in%"uncensored"
                }else{
                    outcome_free <- rep(TRUE,nrow(protocol_data))
                    uncensored <- rep(TRUE,nrow(protocol_data))
                }
                # fit outcome regression
                if (speed & !inherits(try(fit_last <- speedglm::speedglm(formula = interval_outcome_formula,
                                                                         data = protocol_data[outcome_free&uncensored],
                                                                         family = binomial(),
                                                                         maxit = 100),silent = TRUE),
                                      "try-error")){
                } else
                    fit_last <- glm(formula = interval_outcome_formula,
                                    data = protocol_data[outcome_free&uncensored],
                                    family = "binomial")
                # save fitted object
                x$models[["outcome"]][[j]]$fit <- fit_last
                # intervene according to protocol for targets
                # FIXME: intervene on all variables or only those after
                #        time j? those in current outcome_formula
                intervened_data = intervene(
                    ## formula = interval_outcome_formula,
                    ## data = protocol_data[outcome_free&uncensored],
                    data = protocol_data,
                    intervention_table = intervention_table,
                    time = j)
                # set predicted value as outcome for next regression
                ## FIXME: if (target$estimator == "tmle")
                y <- predict(fit_last,newdata = intervened_data ,type = "response")
                current_cnode = protocol_data[[paste0(x$name_censoring,"_",j)]]
                if (length(x$targets[[target_name]]$estimator) == 0 || x$targets[[target_name]]$estimator == "tmle"){
                    W = update_Q(Y = Yhat,
                                 logitQ = lava::logit(y),
                                 cum.g = x$cumulative_intervention_probs[[protocol_name]][,protocol_Cnodes[[j]]], 
                                 uncensored_undeterministic = outcome_free & (current_cnode%in%"uncensored"),
                                 intervention.match = x$intervention_match[[protocol_name]][,intervention_table[time == j-1]$variable])
                }else{
                    W <- predict(fit_last,
                                 newdata = intervened_data ,
                                 type = "response")
                }
                # calculate contribution to influence function
                h.g.ratio <- 1/x$cumulative_intervention_probs[[protocol_name]][,match(paste0("Censored_",j),colnames(x$cumulative_intervention_probs[[protocol_name]]))]
                index <- (current_cnode%in%"uncensored") & intervention_match[,intervention_table[time == j-1]$variable]
                if (any(h.g.ratio[index] != 0)) {
                    x$IC[[target_name]][[protocol_name]][index] <- x$IC[[target_name]][[protocol_name]][index] + (Yhat[index] - W[index]) * h.g.ratio[index]
                }
                ## curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio,
                ## uncensored, intervention.match, regimes.with.positive.weight)
                # prepare next iteration
                set(x = protocol_data,
                    i = which(outcome_free&uncensored), # atrisk: not censored, event-free, no competing risk
                    j = "predicted_outcome",
                    value = W[which(outcome_free&uncensored)])
                # for those who have had an event or died or censored earlier
                set(x = protocol_data,
                    i = which(!outcome_free|!uncensored), # not atrisk
                    j = "predicted_outcome",
                    # fixme j or j-1?
                    value = protocol_data[[paste0(x$name_outcome,"_",j)]][which(!outcome_free|!uncensored)])
            }
            # g-formula and tmle estimator
            x$estimate[[target_name]][[protocol_name]] <- mean(protocol_data$predicted_outcome)
            x$IC[[target_name]][[protocol_name]] <- x$IC[[target_name]][[protocol_name]] + protocol_data$predicted_outcome - mean(protocol_data$predicted_outcome)
        }
    }
    x
}
######################################################################
### run_rtmle.R ends here
