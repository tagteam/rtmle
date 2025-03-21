### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Mar 21 2025 (11:51) 
##           By: Thomas Alexander Gerds
##     Update #: 65
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
##' \item \code{treatment_variables}: A vector with the
##' name(s) of the variable(s) that the protocols intervenes upon. In longitudinal 
##' settings, when a treatment variable is named "A" then the prepared data will contain
##' a column for each of the variables A_0,A_1,----,A_k. This argument can be left unspecified
##' in which case the argument \code{intervention} must be a data.frame with names of treatment variables.
##' See also details.
##' \item \code{intervention}: A vector, or a named data.frame (tibble, data.table), or a function.
##' If it is a vector, it should consist of values 0 and 1 corresponding to what the intervention
##' would set for the \code{treatment_variables}. If it is a data.frame (tibble, data.table),
##' the names specify the treatment variables and the argument \code{treatment_variables} is ignored.
##' Each column must be a factor with levels specifying the treatment options. In longitudinal settings,
##' there can either be only one row in the data.frame and then it is assumed that the intervention is
##' static throughout the followup period. The data.frame can also have a column named time and
##' specify for each time point the values that the intervention sets for the treatment variables.
##' If it is a function, it will be called from \code{intervention_probabilities} with two arguments:
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
            length(value$treatment_variables) == length(intervention_table) &&
            all(intervention_table %in% c(0,1))){
            treatment_variables <- value$treatment_variables
            intervention_table <- data.table::as.data.table(do.call(cbind,lapply(1:length(treatment_variables),function(v){
                factor(intervention_table[[v]],levels = c(0,1))
            })))
            data.table::setnames(intervention_table,treatment_variables)
            x$names$treatment_options <- lapply(treatment_variables,function(x)c(0,1))
            names(x$names$treatment_options) <- treatment_variables
        }else{
            stop("Intervention must be a vector of 0s and 1s or a data.frame (or data.table or tibble) containing factors.")
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
