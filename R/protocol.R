### protocol.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: apr 29 2026 (07:33) 
##           By: Thomas Alexander Gerds
##     Update #: 139
#----------------------------------------------------------------------
##
### Commentary: 
##
### Change Log:
#----------------------------------------------------------------------
## ### Code:
##' Define a treatment protocol for an emulated trial
##'
##' Adds a named protocol to an existing \code{rtmle} object. A protocol
##' defines the values of the treatment variable(s) at each time point
##' during follow-up, including time zero (baseline).
##'
##' @param x An \code{rtmle} object as returned by \code{\link{rtmle_init}}.
##' @param name Name of the protocol.
##' @param intervention A vector, data frame, tibble, or data table specifying
##'   the treatment values dictated by the protocol. If a vector is supplied, it
##'   should contain 0/1 values corresponding to \code{treatment_variables}. If a
##'   data frame, tibble, or data table is supplied, its treatment columns must
##'   be factors whose levels specify the treatment options; in this case
##'   \code{treatment_variables} is ignored. In longitudinal settings, include a
##'   \code{"time"} column in addition to the treatment variables. Ideally,
##'   \code{intervention} supplies values for every value in
##'   \code{x$intervention_nodes}. If fewer rows are supplied and \code{expand =
##'   TRUE}, the last supplied value is used for later time intervals. If
##'   \code{expand = FALSE}, missing time intervals are interpreted as having no
##'   intervention.
##' @param expand Logical. If \code{FALSE} and \code{intervention} contains a
##'   \code{"time"} column, do not expand static interventions across time.
##' @param treatment_variables A vector with the name(s) of the
##'     variable(s) that the protocol intervenes on. In
##'     longitudinal settings, when a treatment variable is named \code{"A"},
##'     the prepared data contain one column for each of
##'     \code{A_0}, \code{A_1}, \dots, \code{A_k}. This argument can be left
##'     unspecified in which case the argument \code{intervention}
##'     must be given.
##' @param intervene_function A character string naming the function used to
##'     intervene under the protocol. Defaults to \code{"intervene"}, which
##'     implements static interventions. The function is called from the internal
##'     function \code{intervention_probabilities()} with two arguments: the
##'     current time interval and the current history of all
##'     variables. It determines the treatment value(s) under the intervention
##'     and should return a matrix with as many columns as there are treatment
##'     variables.
##' @param verbose Logical. If \code{FALSE} suppress all messages. \code{TRUE} is the default.
##' @param ... Not used.
#' @return The modified object contains the treatment variables and
#'     \code{intervention_table} as list elements of \code{x$protocols[[name]]}.
#' @seealso \code{\link{rtmle_init}}, \code{\link{prepare_rtmle_data}},
#'   \code{\link{intervention_match}}, \code{\link{target}},
#'   \code{\link{model_formula}}, \code{\link{run_rtmle}}
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
#' # Intervening on more than one treatment variable
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
            if (length(x$intervention_nodes)>1){
                stop("Argument intervention needs a time variable with values that are equal to or a subset of x$intervention_nodes.")
            }
        }
        if (any(grepl(pattern = "_[0-9]+$",x = treatment_variables))){
            stop("Treatment variables should be given without time suffix.")
        }
        if (any(is.na(too_many <- which(is.na(match(intervention_table[["time"]],x$intervention_nodes,nomatch = NA))))))
            stop(paste0("The following time points are not registered as intervention nodes in the object:",
                        paste0(intervention_table[["time"]][too_many],collapse = ", ")))
        treatment_options <- sapply(treatment_variables,
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
    # check adherence
    x <- intervention_match(x,protocol_name = name)
    x
}

######################################################################
### protocol.R ends here
