### summary.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 29 2024 (10:44) 
## Version: 
## Last-Updated: Sep 22 2024 (18:27) 
##           By: Thomas Alexander Gerds
##     Update #: 17
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Summarizing the results of a register data analysis with the targeted minimum loss estimator 
#'
#' 
#' @param object Object to be summarized
#' @param ...
#' @method summary rtmle
#' @export
summary.rtmle <- function(object,...){
    do.call(rbind,lapply(names(object$targets),function(target_name){
        target <- object$targets[[target_name]]
        protocols <- object$protocols
        do.call(rbind,lapply(names(protocols),function(protocol_name){
            e = object$estimate[[target_name]][[protocol_name]]
            ic = object$IC[[target_name]][[protocol_name]]
            se = sqrt(var(ic)/NROW(ic))
            lower = e-qnorm(.975)*se
            upper = e+qnorm(.975)*se
            sline = Publish::formatCI(x = e,
                                      lower = lower,
                                      upper = upper,
                                      show.x = TRUE)
            data.table(Target = target_name,
                       Protocol = protocol_name,
                       Estimate = e,
                       std.err =se ,
                       "Estimate (CI_95)" = sline)
        }))}))
}


######################################################################
### summary.rtmle.R ends here
