### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Nov  25 2024
##           By: Alessandra
##     Update #: 19
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log: adapt it to have a stochastic intervention
##              In case of stochastic intervention, we need to define the intervene function
#               defined as the probability to receive treatment (A=1) at each time
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
  stopifnot(all(c("name","treatment_variables","intervention","intervention_type")%in% names(value)))
  intervention_times <- x$time[-length(x$time)]
  if( value$intervention_type=="stochastic"){
    stopifnot(length(value$intervene_function)>0)
    value$intervention<- NULL # we do not need to define this at this point
  }
  x$protocols[[value$name]]$intervention_table <- data.table(time = intervention_times,
                                                             variable = paste0(value$treatment_variables,"_",intervention_times),
                                                             value = value$intervention)
  x$protocols[[value$name]]$treatment_variables <- value$treatment_variables
  if (length(value$intervene_function)>0){
    x$protocols[[value$name]]$intervene_function <- value$intervene_function # here we define the function for the stochastic
  }else{
    x$protocols[[value$name]]$intervene_function <- "intervene"
  }
  # we have to save somewhere which kind of intervention we are considering
  x$protocols[[value$name]]$intervention_type<-value$intervention_type
  x
}
######################################################################
### protocol.R ends here








