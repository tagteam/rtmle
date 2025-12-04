### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: dec  3 2025 (06:22) 
##           By: Thomas Alexander Gerds
##     Update #: 94
#----------------------------------------------------------------------
##
### Commentary: 
##
### Change Log:
#----------------------------------------------------------------------
## ### Code:
##' Define a named protocol for an emulated trial
##'
##' This function adds a protocol to an existing object.  A protocol
##' defines the values of the treatment variable(s) at each time point
##' during followup including at time zero (baseline).
##' @param x An rtmle object as obtained with \code{rtmle_init}.
##' @param name The name of the protocol of an emulated trial
##' @param intervention A vector or a data.frame (or tibble or data.table) with names of treatment variables set
##'        to the value that the protocol of an emulated trial dictates. 
##'        If it is a vector, it should consist of values 0 and 1 corresponding to what the
##'     intervention would set for the \code{treatment_variables}. If
##'     it is a data.frame (tibble, data.table), the names specify the
##'     treatment variables and the argument
##'     \code{treatment_variables} is ignored.  Each column must be a
##'     factor with levels specifying the treatment options. In
##'     longitudinal settings, there can either be only one row in the
##'     data.frame and then it is assumed that the intervention is
##'     static throughout the followup period. The data.frame can also
##'     have a column named time and specify for each time point the
##'     values that the intervention sets for the treatment variables.
##' @param treatment_variables A vector with the name(s) of the
##'     variable(s) that the protocols intervenes upon. In
##'     longitudinal settings, when a treatment variable is named "A"
##'     then the prepared data will contain a column for each of the
##'     variables A_0,A_1,----,A_k. This argument can be left
##'     unspecified in which case the argument \code{intervention}
##'     must be given.
##' @param intervene_function A character string: The name of a function
##'        used to intervene under the protocol. Defaults to \code{"intervene"}
##'        which implements static interventions. The function will be called from
##'     the internal function intervention_probabilities with two arguments: the
##'     current time interval and the current history of all
##'     variables. The function determines the value(s) of the
##'     treatment variable(s) under the intervention and should return
##'     a matrix with as many columns as there are treatment variables
##' @param verbose Logical. If \code{FALSE} suppress all messages. \code{TRUE} is the default.
##' @param ... Not (yet) used
#' @return The modified object contains the treatment variables and
#'     the intervention_table as list elemens of
#'     \code{x$protocols[[name]]} where name is given by
#'     \code{value$name}
#' @author Thomas A Gerds \email{tag@@biostat.ku.dk}
#' @examples
#' # ------------------------------------------------------------------------------------------
#' # Intervening on a single treatment variable
#' # ------------------------------------------------------------------------------------------
#' x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- protocol(x,name = "Always_A",
#'                 intervention = data.frame("A" = factor("1",levels = c("0","1"))))
#' x <- protocol(x,name = "Never_A",
#'               intervention = data.frame("A" = factor("0",levels = c("0","1"))))
#' x <- protocol(x,name = "Initiate_A_then_stop",
#'               intervention = data.frame("A" = factor(c("1","0","0"),levels = c("0","1"))))
#' x$protocols
#' # ------------------------------------------------------------------------------------------
#' # Intervening on a more than one treatment variable
#' # ------------------------------------------------------------------------------------------
#' x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- protocol(x,name = "Always_A_never_B",
#'                 intervention = data.frame("A" = factor("1",levels = c("0","1")),
#'                                      "B" = factor("0",levels = c("0","1"))))
#' x <- protocol(x,name = "Always_A_and_B_never_C",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1")),
#'                                               "B" = factor("1",levels = c("0","1")),
#'                                               "C" = factor("0",levels = c("0","1"))))
#' x$protocols
##' @export
protocol <- function(x,
                     name,
                     intervention,
                     treatment_variables,
                     intervene_function = NULL,
                     verbose = TRUE,
                     ...) {
    variable <- NULL
    intervention_times <- x$time[-length(x$time)]
    intervention_table <- intervention
    if (inherits(intervention_table,"data.frame")){
        data.table::setDT(intervention_table)
        treatment_variables <- names(intervention_table)
        if (any(grepl(pattern = "_[0-9]+$",
                      x = treatment_variables))){
            stop("Treatment variables should be given without time suffix.")
        }
        missing_nodes <- (length(intervention_times)-NROW(intervention_table))
        if (missing_nodes>0){
            if (verbose[[1]] == TRUE){
                message("The object specifies more intervention nodes than there are rows in the provided intervention table.\nApply last value carried forward for now, but please check 'x$protocol$intervention_table'.")
            }
            intervention_table <- intervention_table[c(1:NROW(intervention_table),rep(NROW(intervention_table),missing_nodes))]
        }else{
            if (missing_nodes<0){
                if (verbose[[1]] == TRUE){
                    message("The object specifies fewer intervention nodes than there are rows in the provided intervention table.\nCutting these for now, but please check 'x$protocol$intervention_table'.")
                }
                intervention_table <- intervention_table[c(1:length(intervention_times))]
            }
        }
        treatment_options <- sapply(treatment_variables,
                                    function(v){
                                        if (is.factor(intervention[[v]])){
                                            levels(intervention[[v]])
                                        } else{
                                            # FIXME: this seems out of control and will not work
                                            stop("The treatment variables must be factors") 
                                            unique(intervention[[v]])
                                        }
                                    },simplify = FALSE)
    }else{
        if (!missing(treatment_variables) && 
            length(treatment_variables) == length(intervention_table) &&
            all(intervention_table %in% c(0,1))){
            intervention_table <- data.table::as.data.table(lapply(1:length(treatment_variables),function(v){
                factor(intervention_table[[v]],levels = c(0,1))
            }))
            data.table::setnames(intervention_table,treatment_variables)
            treatment_options <- sapply(treatment_variables,function(x)c(0,1),simplify = FALSE)
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
    if (length(intervene_function)>0){
        x$protocols[[name]]$intervene_function <- intervene_function
    }else{
        x$protocols[[name]]$intervene_function <- "intervene"
    }
    x$protocols[[name]]$treatment_variables <- treatment_variables
    x$protocols[[name]]$intervention_table <- intervention_table
    # adding the treatment options if necessary
    if (length(x$names$treatment_options) == 0){
        x$names$treatment_options <- treatment_options
    }else{
        new_options <- setdiff(names(treatment_options),
                               names(x$names$treatment_options))
        if (length(new_options)>0){
            x$names$treatment_options <- c(x$names$treatment_options,treatment_options[new_options])
        }
    }
    x
}

######################################################################
### protocol.R ends here
