### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Apr 25 2025 (10:27)
##           By: Thomas Alexander Gerds
##     Update #: 76
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
##' An optional argument is
##' \itemize{
##' \item \code{verbose}: if FALSE suppress messages
##' }
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
#' @return The modified object contains the treatment variables and the intervention_table
#'         as list elemens of \code{x$protocols[[name]]} where name is given by \code{value$name}
#' @author  Thomas A Gerds \email{tag@@biostat.ku.dk}
#' @examples
#' # ------------------------------------------------------------------------------------------
#' # Intervening on a single treatment variable
#' # ------------------------------------------------------------------------------------------
#' x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' protocol(x) <- list(name = "Always_A",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1"))))
#' protocol(x) <- list(name = "Never_A",
#'                     intervention = data.frame("A" = factor("0",levels = c("0","1"))))
#' protocol(x) <- list(name = "Initiate_A_then_stop",
#'                     intervention = data.frame("A" = factor(c("1","0","0"),levels = c("0","1"))))
#' x$protocols
#' # ------------------------------------------------------------------------------------------
#' # Intervening on a more than one treatment variable
#' # ------------------------------------------------------------------------------------------
#' x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' protocol(x) <- list(name = "Always_A_never_B",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1")),
#'                                               "B" = factor("0",levels = c("0","1"))))
#' protocol(x) <- list(name = "Always_A_and_B_never_C",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1")),
#'                                               "B" = factor("1",levels = c("0","1")),
#'                                               "C" = factor("0",levels = c("0","1"))))
#' x$protocols
##' @export
"protocol<-" <- function(x,...,value) {
  variable <- NULL
  stopifnot(is.list(value))
  stopifnot(all(c("name","intervention") %in% names(value)))
  if (length(value$verbose) == 0) value$verbose <- TRUE
  intervention_times <- x$time[-length(x$time)] # this is defined by the rtmle_initi when you define the intervals
  # in case of a stochastic/dynamic intervention we do not define the intervention
  # as a vector or a data.frame but rather as a function
  # so we will have the definition of the treatment_variables as a factor with 2 levels
  intervention_table <- value$intervention
  if (inherits(intervention_table,"data.frame")){ ## define the intervention as a data.frame with one row for each time for A or A and B
    data.table::setDT(intervention_table)
    treatment_variables <- names(intervention_table)
    if (any(grepl("_[0-9]+$","",treatment_variables))){
      stop("Treatment variables should be given without time suffix.")
    }
    missing_nodes <- (length(intervention_times)-NROW(intervention_table))
    if (missing_nodes>0){
      if (value$verbose[[1]] == TRUE){
        message("The object specifies more intervention nodes than there are rows in the provided intervention table.\nApply last value carried forward for now, 'x$protocol$intervention_table'.")
      }
      intervention_table <- intervention_table[c(1:NROW(intervention_table),rep(NROW(intervention_table),missing_nodes))]
    }else{
      if (missing_nodes<0){
        if (value$verbose[[1]] == TRUE){
          message("The object specifies fewer intervention nodes than there are rows in the provided intervention table.\nCutting these for now, but please check 'x$protocol$intervention_table'.")
        }
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
  }
  else{
    ## this is when you do not define the data.frame but only the treatment_variables name
    ## and the intervention is a vector with one value for each time
    ## this work only with binary treatment
    if ("treatment_variables" %in% names(value) && !is.function(intervention_table) &&
        length(value$treatment_variables) == length(intervention_table) &&
        all(intervention_table %in% c(0,1))){
      treatment_variables <- value$treatment_variables
      # create the intervention table from the vector as a data.table with one raw for each time
      intervention_table <- data.table::as.data.table(lapply(1:length(treatment_variables),function(v){
        factor(intervention_table[[v]],levels = c(0,1))
      }))
      data.table::setnames(intervention_table,treatment_variables)
      x$names$treatment_options <- lapply(treatment_variables,function(x)c(0,1)) ## assigning the two levels myself, but if it is more then two?
      names(x$names$treatment_options) <- treatment_variables
    }else{
      if("treatment_variables" %in% names(value) && is.function(intervention_table))
      { ## we have a stochastic or a dynamic intervention, so it is not fixed, it will depend from the history
        ## and we have to take track of the needed variables for the definition of the intervention
        ## we put something here, but in principle the value is not needed in this case.
        ## to define the levels later on
        treatment_variables<-value$treatment_variables
        # x$names$treatment_options <- lapply(treatment_variables,
        #                                     function(v){treatment_variables[[v]]})
        # names(x$names$treatment_options) <-treatment_variables
        x$names$treatment_options <- lapply(treatment_variables,function(x)c(0,1)) ## assigning the two levels myself, but if it is more then two?
        names(x$names$treatment_options) <- treatment_variables
        intervention_table<-data.table::as.data.table(lapply(1:length(treatment_variables),function(v){
          #factor(NA, levels=c(0,1))
          value$intervention
        }))
        data.table::setnames(intervention_table,treatment_variables)
      }else{stop("Intervention must be a vector of 0s and 1s or a data.frame (or data.table or tibble) containing factors or a function for dynamic or stochastic rules.")}

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

  if (is.function(value$intervention)){
    ## the function inside intervene can be used for a dynamic intervention (intervention value depends on the past values)
    ## or for a stochastic intervention (define a function for the probability of getting the intervention)
    ## these are different because in the dynamic intervention then we get a value for the intervention_table
    ## while for the stochastic intervention you get a probability
    x$protocols[[value$name]]$intervene_function <- value$intervention
    ## and in this case, we should ask for the intervention_type: dynamic or stochastic so to distinguish what to do with it:
    if((value$intervention_type) %in% c("stochastic","dynamic")){
      x$protocols[[value$name]]$intervention_type<-value$intervention_type
    }
    else{
      stop("you defined a function as for the intervention, you need to specify the intervention_type as \"stochastic \" or \" dynamic \" ")
    }
    ## in case of a function we should define the variables the rule depends on
    ## DISCUSS WITH THOMAS AND UNDERSTAND HOW TO DO THIS, AS HE DID FOR THE gform and Qform (target)
  }else{
    x$protocols[[value$name]]$intervene_function <- "intervene"
    x$protocols[[value$name]]$intervention_type<-"static"
  }
  x$protocols[[value$name]]$treatment_variables <- treatment_variables
  x$protocols[[value$name]]$intervention_table <- intervention_table
  x
}
