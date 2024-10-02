### coef.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:48) 
## Version: 
## Last-Updated: Oct  2 2024 (09:32) 
##           By: Thomas Alexander Gerds
##     Update #: 5
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
coef.rtmle <- function(object,...){
    res = lapply(names(object$models),function(p){
        inner = lapply(names(object$models[[p]]),function(m){
        process = lapply(names(object$models[[p]][[m]]),function(j){
            if (length(object$models[[p]][[m]][[j]]$fit)>0)
                ## coef(summary(object$models[[p]][[m]][[j]]$fit))
                Publish::publish(object$models[[p]][[m]][[j]]$fit)
            else
                NULL
        })
        names(process) = names(object$models[[p]][[m]])
        process
        })
    })
    names(res) = names(object$models)
    res
}
######################################################################
### coef.rtmle.R ends here
