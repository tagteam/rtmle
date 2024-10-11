### summary.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 29 2024 (10:44) 
## Version: 
## Last-Updated: Oct 11 2024 (08:15) 
##           By: Thomas Alexander Gerds
##     Update #: 36
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
#' @param ... not used
#' @export
#' @method summary rtmle
summary.rtmle <- function(object,...){
    Estimate <- Upper <- Lower <- NULL
    do.call(rbind,lapply(names(object$targets),function(target_name){
        target <- object$targets[[target_name]]
        protocols <- object$protocols
        out <- do.call(rbind,lapply(names(protocols),function(protocol_name){
            # to avoid the internal selfdetect problem we take a copy
            e = data.table::copy(object$estimate[[target_name]][[protocol_name]])
            e[, "Estimate (CI_95)":= Publish::formatCI(x = 100*Estimate,lower = 100*Lower,upper = 100*Upper,show.x = TRUE)]
        }))
        for (nix in c("Estimate","Standard_error","Lower","Upper"))
            set(out,j = nix,value = NULL)
        out
    }))
}


######################################################################
### summary.rtmle.R ends here
