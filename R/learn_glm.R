### rtmle_glm.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (12:49) 
## Version: 
## Last-Updated: Sep 26 2024 (09:37) 
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
learn_glm <- function(formula,data, speed = TRUE,...){
    speed = FALSE
    if (speed && !inherits(try(
                      fit <- speedglm::speedglm(formula = formula,data = data,family = binomial(),maxit = 100),silent = TRUE),
                      "try-error")){
    } else{
        if (inherits(try(fit <- glm(formula = formula,data = data,family = binomial()),
                         silent = TRUE),"try-error")){
            ff = as.character(formula)
            stop(paste0("Could not fit this model with glm:\n","Outcome: ",ff[[2]],"\nRight hand side: ",ff[[3]]))
        }
    }
    class(fit) <- c(class(fit),"glm")
    fit
}


######################################################################
### rtmle_glm.R ends here
