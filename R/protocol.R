### protocol.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: mar 27 2026 (06:25) 
##           By: Thomas Alexander Gerds
##     Update #: 135
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
##'     it is a data.frame (tibble, data.table), the names of the data.frame are the
##'     treatment variables. (In this case the argument
##'     \code{treatment_variables} is ignored.) Each column of the data.frame must be a
##'     factor with levels specifying the treatment options. In
##'     longitudinal settings, the data.frame should have column called \code{"time"} in
##'     addition to the treatment variables. Ideally the data.frame \code{intervention} provides values
##'     for all times stored in \code{x$intervention_nodes}. However, there are two special cases.
##'     If the number of rows in the data.frame \code{intervention} is smaller than
##'    \code{length(x$intervention_nodes)}, then 1. if argument \code{expand} is \code{TRUE},
##'    the last value (the only value if there is only one row) is set for all later time intervals.
##'     2. if argument \code{expand} is \code{FALSE} it is assumed that no intervention is wanted at
##'    the time intervals that are missing in data.frame \code{intervention}.
##' @param expand Logical. If \code{FALSE} and time is a column in the data.frame
##'     given by argument \code{intervention}, then do not expand static interventions across time. 
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
#' x <- rtmle_init(time_grid=0:3,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- protocol(x,name = "Always_A",
#'                 intervention = data.frame(time=x$intervention_nodes,
#'                                           "A" = factor("1",levels = c("0","1"))))
#' x <- protocol(x,name = "Never_A",
#'               intervention = data.frame(time=x$intervention_nodes,
#'                                         "A" = factor("0",levels = c("0","1"))))
#' x <- protocol(x,name = "Initiate_A_then_stop",
#'               intervention = data.frame(time=x$intervention_nodes,
#'                                         "A" = factor(c("1","0","0"),levels = c("0","1"))))
#' x$protocols
#' # ------------------------------------------------------------------------------------------
#' # Intervening on a more than one treatment variable
#' # ------------------------------------------------------------------------------------------
#' x <- rtmle_init(time_grid=0:3,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- protocol(x,name = "Always_A_never_B",
#'                 intervention = data.frame(time=x$intervention_nodes,
#'                                      "A" = factor("1",levels = c("0","1")),
#'                                      "B" = factor("0",levels = c("0","1"))))
#' x <- protocol(x,name = "Always_A_and_B_never_C",
#'                     intervention = data.frame(time=x$intervention_nodes,
#'                                               "A" = factor("1",levels = c("0","1")),
#'                                               "B" = factor("1",levels = c("0","1")),
#'                                               "C" = factor("0",levels = c("0","1"))))
#' x$protocols
##' @export
protocol <- function(x,
                     name,
                     intervention,
                     expand = TRUE,
                     treatment_variables,
                     intervene_function = NULL,
                     verbose = TRUE,
                     ...) {
    variable <- NULL
    #
    # User provided a data.frame
    #
    if (inherits(intervention,"data.frame")){
        intervention_table <- data.table::copy(intervention)
        data.table::setDT(intervention_table)
        treatment_variables <- names(intervention_table)
        if ((time_position <- match("time",treatment_variables,nomatch = 0))>0){
            treatment_variables <- treatment_variables[-time_position]
        }else{
            if (length(x$time)>1){
                stop("Argument intervention needs a time variable with values that are equal to or a subset of x$intervention_nodes.")
            }
        }
        if (any(grepl(pattern = "_[0-9]+$",x = treatment_variables))){
            stop("Treatment variables should be given without time suffix.")
        }
        if (any(is.na(too_many <- which(is.na(match(intervention_table[["time"]],x$intervention_nodes,nomatch = NA))))))
            stop(paste0("The following time points are not registered as intervention nodes in the object:",
                        paste0(intervention_table[["time"]][too_many],collapse = ", ")))
        # FIXME: do something with missing nodes? we could add NA but not having them may have the same effect
        ## missing_nodes <- match(x$intervention_nodes,intervention_table[["time"]],nomatch = NA)
        treatment_options <-
            sapply(treatment_variables,
                   function(v){
                       if (is.factor(intervention[[v]])){
                           if (length(levels(intervention[[v]])) == 2){
                               levels(intervention[[v]])
                           }else{
                               stop(paste0("All treatment variables must have exactly 2 levels. Problem with variable ",
                                           v),".")
                           }
                       } else{
                           stop(paste0("The treatment variables must be factors. Problem with variable ",v,"."))
                       }
                   },simplify = FALSE)
    }else{
        #
        # User provided a intervention values and separately names of treatment variables 
        # 
        if (!missing(treatment_variables) && 
            length(treatment_variables) == length(intervention) &&
            all(intervention %in% c(0,1))){
            intervention_table <- data.table::as.data.table(lapply(1:length(treatment_variables),function(v){
                factor(intervention[[v]],levels = c(0,1))
            }))
            intervention_table <- cbind(time = x$intervention_nodes,
                                        intervention_table)
            data.table::setnames(intervention_table,c("time",treatment_variables))
            treatment_options <- sapply(treatment_variables,function(x)c(0,1),simplify = FALSE)
        }else{
            stop("Argument `intervention' is not a data.frame. Hence it must be a vector of 0s and 1s
                  with the same length as argument treatment_variables")
        }
    }
    # turn from wide into long format
    intervention_table <- data.table::melt(
                                          intervention_table, 
                                          id.vars = "time", 
                                          variable.name = "variable", 
                                          value.name = "value",
                                          value.factor = TRUE
                                      )[, variable := paste0(variable, "_", time)]
    if (length(intervene_function)>0){
        x$protocols[[name]]$intervene_function <- intervene_function
    }else{
        x$protocols[[name]]$intervene_function <- "intervene"
    }
    x$protocols[[name]]$treatment_variables <- treatment_variables
    x$protocols[[name]]$intervention_table <- intervention_table[]
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
