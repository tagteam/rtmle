### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Dec  6 2024 (13:59) 
##           By: Thomas Alexander Gerds
##     Update #: 26
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
##' \item \code{treatment_variables}: the name(s) of the variable(s) that the protocols intervenes upon
##' \item \code{intervention}: A matrix or a function. If it is a matrix it should contain the values
##' that the variables are set to under the intervention as columns and one row for each time point.
##' If the intervened values are the same for all time points it is sufficient to provide a single value per variable.
##' See examples. If it is a function, it will be called from \code{intervention_probabilities} with two arguments: 
##' the current time interval and the current history of all variables. The function determines the value(s)
##' of the treatment variable(s) under the intervention and should return a matrix with as many columns as there are
##' treatment variables.
##' }
##' @export
"protocol<-" <- function(x,...,value) {
    stopifnot(is.list(value))
    stopifnot(all(c("name","treatment_variables","intervention")%in%names(value)))
    intervention_times <- x$time[-length(x$time)]
    tv <- value$treatment_variables
    if (length(grep("_[0-9]+$",tv))>0){
        varnames <- unique(sub("_[0-9]+$",tv))
    }else{
        varnames <- value$treatment_variables
    }
    if (length(varnames)>1){
        # multiple treatment nodes
        if (is.list(tv)) {
            stopifnot(length(unique(sapply(tv,length))))
            tv <- do.call(cbind,tv)
        }else{
            if (!is.matrix(tv)) 
                tv <- do.call(cbind,lapply(varnames,function(v)paste0(v,"_",intervention_times)))
        }
    }
    else{
        if (lengths(tv) == 1) {
            tv <- paste0(tv,"_",intervention_times)
        } else {
            if (lengths(tv) != lengths(intervention_times)){
                stop(paste0("Argument treatment_variables has the wrong length. You can either provide",
                            "the name of the treatment variable such as 'A' as a character string without the '_k' subscript,",
                            "or the vector for all time points such as 'c('A_0', 'A_1', ..., 'A_k') where k=",length(intervention_times),
                            " is the number of time points including 0 but minus the last time point."))
            }
        }
    }
    it <- cbind(data.table(time = intervention_times),
                tv,
                data.table(value = value$intervention))
    if (length(value$intervene_function)>0){
        x$protocols[[value$name]]$intervene_function <- value$intervene_function
    }else{
        x$protocols[[value$name]]$intervene_function <- "intervene"
    }
    x$protocols[[value$name]]$treatment_variables <- tv
    x$protocols[[value$name]]$intervention_table <- it
    x
}
######################################################################
### protocol.R ends here
