### summary.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 29 2024 (10:44) 
## Version: 
## Last-Updated: Aug 26 2024 (10:07) 
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
            e = target$estimates
            ic = target$IC
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
