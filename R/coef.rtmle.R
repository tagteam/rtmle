### coef.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:48) 
## Version: 
## Last-Updated: Jul 19 2024 (10:50) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
coef.rtmle <- function(object,...){
    res = lapply(names(object$models),function(m){
        process = lapply(names(object$models[[m]]),function(j){
            if (length(object$models[[m]][[j]]$fit)>0)
                coef(summary(object$models[[m]][[j]]$fit))
            else
                NULL
        })
        names(process) = names(object$models[[m]])
        process
    })
    names(res) = names(object$models)
    res
}
######################################################################
### coef.rtmle.R ends here
