### fit_propensity_models.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 29 2024 (14:01) 
## Version: 
## Last-Updated: Jul 29 2024 (14:51) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

fit_propensity_models <- function(x,...){
    for (pro in x$protocols){
        pro <- x$protocols[[x$targets[[tar]]$protocol]]$protocol
    for (j in 1:length(x$models[["propensity"]])){
        # data at-risk at the beginning of the interval
        if (j > 1){
            outcome_free <- x$prepared_data[[paste0(x$name_outcome,"_",j-1)]]%in%0
            # competing risks
            if (length(x$Dnodes)>0) outcome_free <- outcome_free&x$prepared_data[[paste0(x$name_competing,"_",j-1)]]%in%0
            uncensored <- x$prepared_data[[paste0(x$name_propensity,"_",j-1)]]%in%"uncensored"
        }else{
            outcome_free <- rep(TRUE,nrow(x$prepared_data))
            uncensored <- rep(TRUE,nrow(x$prepared_data))
        }
        if (any(outcome_free&uncensored)){
            # fit the propensity regression model
            if (!is.null(ff <- x$models[["propensity"]][[j]]$formula)){
                if (speed & !inherits(try(
                                 x$models[["propensity"]][[j]]$fit <- speedglm::speedglm(formula = ff,data = x$prepared_data[outcome_free&uncensored],family = binomial(),maxit = 100),silent = TRUE),
                                 "try-error")){
                } else{
                    if (inherits(try(x$models[["propensity"]][[j]]$fit <- glm(formula = ff,data = x$prepared_data[outcome_free&uncensored],family = binomial()),
                                     silent = TRUE),"try-error"))
                        stop(paste0("Could not fit ",m," model"))
                }

            }
        }
    }
}

######################################################################
### fit_propensity_models.R ends here
