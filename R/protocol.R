### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Jul 26 2024 (11:57) 
##           By: Thomas Alexander Gerds
##     Update #: 14
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
    x$protocols[[value$name]]$protocol <- data.table(time = x$time,
                                                     variable = paste0(value$treatment_variables,"_",x$time),
                                                     value = value$intervention)
    x$protocols[[value$name]]$treatment_variables <- value$treatment_variables
    x
}
######################################################################
### protocol.R ends here
