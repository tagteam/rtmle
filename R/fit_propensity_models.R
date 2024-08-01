### fit_propensity_models.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 29 2024 (14:01) 
## Version: 
## Last-Updated: Aug  1 2024 (10:45) 
##           By: Thomas Alexander Gerds
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

fit_propensity_models <- function(formula,data, speed = TRUE,...){
    if (speed & !inherits(try(
                     fit <- speedglm::speedglm(formula = formula,data = data,family = binomial(),maxit = 100),silent = TRUE),
                     "try-error")){
    } else{
        if (inherits(try(fit <- glm(formula = formula,data = data,family = binomial()),
                         silent = TRUE),"try-error"))
            stop(paste0("Could not fit propensity score model for protocol ",protocol_name,"."))
    }
    fit
}
######################################################################
### fit_propensity_models.R ends here
