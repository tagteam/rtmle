### model_rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: Oct  2 2024 (15:35) 
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
#' @export
"model<-" <- function(x,...,value) {
    stopifnot(is.list(value))
    stopifnot(all(c("formalizer","treatment_variables")%in%names(value)))
    if (is.character(value$formalizer) && value$formalizer == "additive"){
        x$models = additive_formalizer(x = x,
                                       treatment_variables = value$treatment_variables,
                                       Markov = NULL)
    }
    x
}
######################################################################
### model_rtmle.R ends here
