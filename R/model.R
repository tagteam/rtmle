### model_rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: Jul 26 2024 (11:57) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @export
"model<-" <- function(x,...,value) {
    stopifnot(is.list(value))
    stopifnot(all(c("formalizer","treatment_variables")%in%names(value)))
    if (is.character(value$formalizer) && value$formalizer == "additive"){
        x$models = additive_formalizer(x = x,
                                       treatment_variables = value$treatment_variables,
                                       Markov = NULL)
    }else{
        x$models[[model]][[paste0("time_",time)]]$formula <- value
    }
    x
}
######################################################################
### model_rtmle.R ends here
