### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Dec  4 2024 (14:15) 
##           By: Thomas Alexander Gerds
##     Update #: 20
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
## ### Code:
##' Define a named protocol for a hypothetical/emulated trial
##'
##' This function adds a protocol to an existing object.
##' A protocol defines the values of the treatment variable(s)
##' at each time point during followup including at time zero (baseline).
##' @param x An rtmle object as obtained with \code{rtmle_init}.
##' @param ... Not (yet) used
##' @param value list with three forced elements:
##' \itemize{
##' \item \code{name}: the name of the protocol
##' \item \code{variable}: the name(s) of the variable(s) that the protocols intervenes upon
##' \item \code{intervention}: function which determines the value(s) of the variable(s) under the intervention
##' }
##' @export
"protocol<-" <- function(x,...,value) {
    stopifnot(is.list(value))
    stopifnot(all(c("name","treatment_variables","intervention")%in%names(value)))
    intervention_times <- x$time[-length(x$time)]
    tv <- value$treatment_variables
    if (length(grep("_[0-9]+$",tv))>0){
        varnames <- unique(sub("_[0-9]+$",tv))
        if (length(varnames)>1)
            if (is.list(tv)) {
                stopifnot(length(unique(sapply(tv,length))))
                tv <- do.call(cbind,tv)
            }else{
                tv <- do.call(cbind,lapply(varnames,function(v)paste0(v,"_",intervention_times)))
            }
        it <- cbind(data.table(time = intervention_times),
                    tv,
                    data.table(value = value$intervention))
    }else{
        x$protocols[[value$name]]$treatment_variables <- tv
        x$protocols[[value$name]]$intervention_table <- data.table(time = intervention_times,
                                                                   variable = paste0(value$treatment_variables,"_",intervention_times),
                                                                   value = value$intervention)
        if (length(value$intervene_function)>0){
            x$protocols[[value$name]]$intervene_function <- value$intervene_function
        }else{
            x$protocols[[value$name]]$intervene_function <- "intervene"
        }
    }
    x
}
######################################################################
### protocol.R ends here
