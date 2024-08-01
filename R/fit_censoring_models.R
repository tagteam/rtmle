### fit_censoring_models.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 29 2024 (13:59) 
## Version: 
## Last-Updated: Aug  1 2024 (10:45) 
##           By: Thomas Alexander Gerds
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
fit_censoring_models <- function(formula,data, speed = TRUE,...){
    if (speed & !inherits(try(
                     fit <- speedglm::speedglm(formula = formula,data = data,family = binomial(),maxit = 100),silent = TRUE),
                     "try-error")){
    } else{
        if (inherits(try(fit <- glm(formula = formula,data = data,family = binomial()),
                         silent = TRUE),"try-error"))
            stop(paste0("Could not fit censoring model."))
    }
    fit
}


######################################################################
### fit_censoring_models.R ends here
