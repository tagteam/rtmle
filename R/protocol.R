### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Feb 18 25
##           By: Alessandra
##     Update #: 29
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
##' \item \code{treatment_variables}: A list/Matrix/vector with the name(s) of the variable(s) that the protocols intervenes upon.
##' If only one treatment, then it would be one value with the name of treatment "A" or the vector with one name for each time
##' A_0,A_1,----,A_k.  If more than one treatments, then it can be defined as a list of length equal to the number of treatment
##' where each argument of the list has the vector of the treatment at different times $A: A_0,A_1,...,A_k, $B:B_0,B_1,...B_k
##' OR it can be a matrix where each column contains the variables names of treatment at different times.
##' \item \code{intervention}: A list/matrix/vector or a function. If it is a matrix it should contain the values
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
    stopifnot(all(c("name","treatment_variables","intervention") %in% names(value)))
    intervention_times <- x$time[-length(x$time)]
    tv<-value$treatment_variables

    if(is.list(tv)){
      varnames <- names(tv) # names of the treatments
      stopifnot(length(unique(sapply(tv,length)))==1)
      tv_seq <- do.call(cbind,tv)
    }
    else{
      if(is.matrix(tv)){

        if( dim(tv)[1]>1 & dim(tv)[1] != lengths(intervention_times)){
          stop(paste0("when using a matrix for the definition of the treatment_variables",
                      "then each column has to corespond to one treatment",
                      "and then the number of row needs to be of",length(intervention_times),
                      ",whic is the number of time points including 0 but minus the last time point. "))}
        else{
        tv_seq<-tv
        varnames <- unique(sub("_[0-9]+$","",tv))
        }
      }

      else{ #if it is not a list, nor a matrix then we have only one treatment and
        # check if they provide as input each time OR only 1 value
        if (length(grep("_[0-9]+$",tv))>0){

          if (length(grep("_[0-9]+$",tv)) != lengths(intervention_times)){
            stop(paste0("Argument treatment_variables has the wrong length. You can either provide",
                        "the name of the treatment variable such as 'A' as a character string without the '_k' subscript,",
                        "or the vector for all time points such as 'c('A_0', 'A_1', ..., 'A_k') where k=",length(intervention_times),
                        " is the number of time points including 0 but minus the last time point."))}
          tv_seq<-tv
          varnames <- unique(sub("_[0-9]+$","",tv))
        }
        else{

          if(lengths(tv) != 1){
            stop(paste0("Argument treatment_variables has the wrong length. You can either provide",
                        "the name of the treatment variable such as 'A' as a character string without the '_k' subscript,",
                        "or the vector for all time points such as 'c('A_0', 'A_1', ..., 'A_k') where k=",length(intervention_times),
                        " is the number of time points including 0 but minus the last time point."))}

          varnames <- value$treatment_variables
          tv_seq <- as.vector(sapply(varnames, FUN=function(s){paste0(s,"_",intervention_times)}))
        }

      }



    }




    #---- we look at the intervention now:
    if(is.list(value$intervention)){
      intervention <- do.call(cbind,value$intervention)
      # if we work on the combination--> makes sense to actually provide
      # the values for A and B for the intervention table
      # since we will use this for the predict and they might have two separate effects
    }

    else{

      if(is.matrix(value$intervention)){
        # make sure it is the correct dimension:
        if( dim(tv)[value$intervention]>1 & dim(value$intervention)[1] != lengths(intervention_times)){
          stop(paste0("when using a matrix for the definition of the intervention",
                      "then each column has to corespond to one treatment",
                      "and then the number of row needs to be of",length(intervention_times),
                      ",whic is the number of time points including 0 but minus the last time point. "))}

        #if it has the correct dimension--> save it like that
        intervention <- value$intervention

      }

      ## if it is not a list, nor a matrix then there is only one treatment and it is a number or a vector
      if(length(varnames)>1 | (length(varnames)==1 & length(value$intervention) != length(intervention_times) ) ){
        stop(paste0("Argument intervention has the wrong length.
                    You are either specifying the intervention for only one treatment,
                    while in treatment_variables there are more treatments, OR the length of the vector for the intervention is wrong."))}

      intervention<-ifelse(length(value$intervention)==1, rep(value$intervention, length(intervention_times)), value$intervention)

    }


    it <- cbind(data.table(time = intervention_times),
                variable = tv_seq,
                data.table(value = intervention ))

    if (length(value$intervene_function)>0){
        x$protocols[[value$name]]$intervene_function <- value$intervene_function
    }else{
        x$protocols[[value$name]]$intervene_function <- "intervene"
    }
    x$protocols[[value$name]]$treatment_variables <- varnames
    x$protocols[[value$name]]$intervention_table <- it
    x
}
######################################################################
### protocol.R ends here
