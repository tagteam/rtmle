### run_rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  1 2024 (09:11) 
## Version: 
## Last-Updated: Jul 29 2024 (14:51) 
##           By: Thomas Alexander Gerds
##     Update #: 137
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @export
run_rtmle <- function(x,speed = FALSE,...){
    #
    # g-part: fit nuisance parameter models for propensity and censoring
    #
    # FIXME: for target in targets
    target = 1
    x$targets[[target]]$IC <- numeric(nrow(x$prepared_data))
    intervention_probs <- data.table(ID = x$prepared_data[[x$name_id]])
    intervention_probs <- cbind(intervention_probs,matrix(1,nrow = NROW(x$prepared_data),ncol = length(x$Anodes)+length(x$Cnodes)))
    setnames(intervention_probs,c(x$name_id,c(rbind(x$Anodes,x$Cnodes))))

                # predict the propensity score/1-probability of censored
                # intervene according to protocol for targets
                # FIXME: predict also those censored?
                intervened_data = intervene(formula = ff,
                                            data = x$prepared_data[outcome_free],
                                            protocol = x$protocols[[x$targets[[1]]$protocol]]$protocol,
                                            time = j)
                if (m == "propensity")
                    current_node <- x$Anodes[[j]]
                else
                    current_node <- x$Cnodes[[j]]
                intervention_probs[outcome_free][[current_node]] <- predict(x$models[[m]][[j]]$fit,newdata = intervened_data,type = "response")
    
    x$cumulative_intervention_probs <- matrixStats::rowCumprods(as.matrix(intervention_probs[,-1,with = FALSE]))
    # fixme for target in targets
    target <- 1
    pro <- x$protocols[[x$targets[[target]]$protocol]]$protocol
    intervention_match <- matrix(0,
                                 ncol = length(x$Anodes),
                                 nrow = nrow(x$prepared_data))
    for(j in x$times[-c(length(x$times))]){
        if (j == 0)
            intervention_match[,j+1] <- previous <- (x$prepared_data[[pro[j+1]$variable]] %in% c(pro[j+1]$value,NA))
        else
            intervention_match[,j+1] <- previous <- previous*(x$prepared_data[[pro[j+1]$variable]] %in% c(pro[j+1]$value,NA))
    }
    colnames(intervention_match) <- x$Anodes
    x$intervention_match = intervention_match
    # 
    # Q-part: loop backwards in time through iterative condtional expectations
    #
    # the first step is the outcome regression of the last time interval
    for (j in rev(x$times)[-length(x$times)]){
        # formula
        interval_outcome_formula = formula(x$models[["outcome"]][[paste0(x$name_outcome,"_",j)]]$formula)
        #  at the last time interval the observed outcome is used
        if (j != max(x$times)) {
            interval_outcome_formula <- update(interval_outcome_formula,"predicted_outcome~.")
            Yhat <- x$prepared_data[["predicted_outcome"]]
        }else{
            Yhat <- x$prepared_data[[all.vars(interval_outcome_formula)[[1]]]]
        }
        ## print(interval_outcome_formula)
        ## print(head(Yhat))
        # data at-risk at the beginning of the interval
        if (j > 1){
            outcome_free <- x$prepared_data[[paste0(x$name_outcome,"_",j-1)]]%in%0
            # competing risks
            if (length(x$Dnodes)>0) outcome_free <- outcome_free&x$prepared_data[[paste0(x$name_competing,"_",j-1)]]%in%0
            uncensored <- x$prepared_data[[paste0(x$name_censoring,"_",j-1)]]%in%"uncensored"
        }else{
            outcome_free <- rep(TRUE,nrow(x$prepared_data))
            uncensored <- rep(TRUE,nrow(x$prepared_data))
        }
        # fit outcome regression
        if (speed & !inherits(try(fit_last <- speedglm::speedglm(formula = interval_outcome_formula,
                                                                 data = x$prepared_data[outcome_free&uncensored],
                                                                 family = binomial(),
                                                                 maxit = 100),silent = TRUE),
                              "try-error")){
        } else
            fit_last <- glm(formula = interval_outcome_formula,
                            data = x$prepared_data[outcome_free&uncensored],
                            family = "binomial")
        # save fitted object
        x$models[["outcome"]][[j]]$fit <- fit_last
        # intervene according to protocol for targets
        intervened_data = intervene(formula = interval_outcome_formula,
                                    ## data = x$prepared_data[outcome_free&uncensored],
                                    data = x$prepared_data,
                                    protocol = pro,
                                    time = j)
        # set predicted value as outcome for next regression
        ## FIXME: if (target$estimator == "tmle")
        y <- predict(fit_last,newdata = intervened_data ,type = "response")
        current_cnode = x$prepared_data[[paste0(x$name_censoring,"_",j)]]
        if (length(x$targets[[target]]$estimator) == 0 || x$targets[[target]]$estimator == "tmle"){
            W = update_Q(Y = Yhat,
                         logitQ = lava::logit(y),
                         cum.g = x$cumulative_intervention_probs[,match(paste0("Censored_",j),colnames(x$cumulative_intervention_probs))], 
                         uncensored_undeterministic = outcome_free & (current_cnode%in%"uncensored"),
                         intervention.match = intervention_match[,pro[time == j-1]$variable])
        }else{
            W <- predict(fit_last,
                         newdata = intervened_data ,
                         type = "response")
        }
        # calculate contribution to influence function
        h.g.ratio <- 1/x$cumulative_intervention_probs[,match(paste0("Censored_",j),colnames(x$cumulative_intervention_probs))]
        index <- (current_cnode%in%"uncensored") & intervention_match[,pro[time == j-1]$variable]
        if (any(h.g.ratio[index] != 0)) {
            x$targets[[target]]$IC[index] <- x$targets[[target]]$IC[index] + (Yhat[index] - W[index]) * h.g.ratio[index]
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
            i = which(!outcome_free|!uncensored), # not atrisk
            j = "predicted_outcome",
            # fixme j or j-1?
            value = x$prepared_data[[paste0(x$name_outcome,"_",j)]][which(!outcome_free|!uncensored)])
    }
    # g-formula and tmle estimator
    x$targets[[target]]$estimate <- mean(x$prepared_data$predicted_outcome)
    x$targets[[target]]$IC <- x$targets[[target]]$IC + x$prepared_data$predicted_outcome - mean(x$prepared_data$predicted_outcome)
    x
}
######################################################################
### run_rtmle.R ends here
