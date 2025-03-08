### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Mar  8 2025 (08:10) 
##           By: Thomas Alexander Gerds
##     Update #: 57
#----------------------------------------------------------------------
##
### Commentary: We want to change this for 2 treatments, where the intervention is defined on both
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
##' \item \code{treatment_variable}: A list/Matrix/vector with the
##' name(s) of the variable(s) that the protocols intervenes upon.
##' If only one treatment, then it would be one value with the name of treatment
##' "A" or the vector with one name for each time
##' A_0,A_1,----,A_k.  If more than one treatments, then it can be defined as a list
##' of length equal to the number of treatment
##' where each argument of the list has the vector of the treatment at different times
##' $A: A_0,A_1,...,A_k, $B:B_0,B_1,...B_k
##' OR it can be a matrix where each column contains the variables names
##' of treatment at different times.
##' \item \code{intervention}: A list/matrix/vector or a function.
##'
##' If it is a matrix it should contain the values
##' that the variables are set to under the intervention as columns
##' and one row for each time point.
##' If the intervened values are the same for all time points it is sufficient to provide a single value per variable.
##' See examples. If it is a function, it will be called from \code{intervention_probabilities} with two arguments:
##' the current time interval and the current history of all variables. The function determines the value(s)
##' of the treatment variable(s) under the intervention and should return a matrix with as many columns as there are
##' treatment variables.
##' }
##' @export
"protocol<-" <- function(x,...,value) {
    variable <- NULL
    stopifnot(is.list(value))
    stopifnot(all(c("name","intervention") %in% names(value)))
    intervention_times <- x$time[-length(x$time)]
    intervention_table <- value$intervention
    if (inherits(intervention_table,"data.frame")){
        data.table::setDT(intervention_table)
        treatment_variables <- names(intervention_table)
        if (any(grepl("_[0-9]+$","",treatment_variables))){
            stop("Treatment variables should be given without time suffix.")
        }
        missing_nodes <- (length(intervention_times)-NROW(intervention_table))
        if (missing_nodes>0){
            warning("The object specifies more intervention nodes than there are rows in the provided intervention table.\nApply last value carried forward for now.")
            intervention_table <- intervention_table[c(1:NROW(intervention_table),rep(NROW(intervention_table),missing_nodes))]
        }else{
            if (missing_nodes<0){
                warning("The object specifies fewer intervention nodes than there are rows in the provided intervention table.\nCutting these for now.")
                intervention_table <- intervention_table[c(1:length(intervention_times))]
            }
        }
        x$names$treatment_options <- lapply(treatment_variables,
                                            function(v){
                                                if (is.factor(value$intervention[[v]])){
                                                    levels(value$intervention[[v]])
                                                } else{
                                                    # FIXME: this seems out of control and will not work
                                                    stop("The treatment variables must be factors") 
                                                    unique(value$intervention[[v]])
                                                }
                                            })
        names(x$names$treatment_options) <- treatment_variables
    }else{
        if ("treatment_variables" %in% names(value) &&
            length(value$treatment_variables) == 1  &&
            length(intervention_table) == 1 &&
            intervention_table %in% c(0,1)){
            intervention_table <- data.table::as.data.table(factor(intervention_table,levels = c(0,1)))
            names(intervention_table) <- value$treatment_variables
            treatment_variables <- value$treatment_variables
            x$names$treatment_options <- list(c(0,1))
            names(x$names$treatment_options) <- treatment_variables
        }else{
            stop("Intervention must be a data.frame (or data.table or tibble)")
        }
    }
    if (!("time" %in% names(intervention_table))){
        intervention_table <- cbind(time = intervention_times,
                                    intervention_table)
    }else{
        if (!(all(intervention_table[["time"]] == intervention_times))){
            stop(paste("Intervention times do not match. Object contains:\n",
                       paste(intervention_times,collapse = ", ")))
        }
    }
    intervention_table <- do.call("rbind",lapply(treatment_variables,function(v){
        itv <- intervention_table[,c("time",v),with = FALSE]
        itv[,variable := paste0(v,"_",itv$time)]
        data.table::setnames(itv,old = v,new = "value")
        data.table::setcolorder(itv,c("time","variable","value"))
        itv
    }))
    if (length(value$intervene_function)>0){
        x$protocols[[value$name]]$intervene_function <- value$intervene_function
    }else{
        x$protocols[[value$name]]$intervene_function <- "intervene"
    }
    x$protocols[[value$name]]$treatment_variables <- treatment_variables
    x$protocols[[value$name]]$intervention_table <- intervention_table
    x
}
######################################################################
### protocol.R ends here
