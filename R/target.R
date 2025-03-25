### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Dec 12 2024 (18:34)
##           By: Thomas Alexander Gerds
##     Update #: 32
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
## ### Code:
##' Define a target (parameter) of a hypothetical/emulated trial
##'
##' A target defines a target parameter, an estimator and models
##' for the nuisance parameters.
##' @param x An rtmle object as obtained with \code{rtmle_init}.
##' @param ... Not (yet) used
##' @param value list with three forced elements:
##' \itemize{
##' \item \code{name}: the name of the target parameter.
##' \item \code{strategy}: the nuisance parameter modeling strategy consists of formulas and libraries for the regression models.
##' \item \code{estimator}: the estimator of the target parameter.
##' }
##' @export
"target<-" <- function(x,...,value) {
    stopifnot(is.list(value))
    stopifnot(all(c("name","strategy","estimator","protocols")%in%names(value)))
    stopifnot(all(value[["protocols"]]%in%names(x$protocols)))
    if (length(value[["include_variables"]])>0)
        include_variables <- value[["include_variables"]]
    else
        include_variables <- NULL
    if (length(value[["exclude_variables"]])>0)
        exclude_variables <- value[["exclude_variables"]]
    else
        exclude_variables <- NULL
    if (!is.null(x$continuous_outcome)){
      continuous_outcome <- x$continuous_outcome
    } else{
      continuous_outcome <- FALSE
    }
    all_treatment_variables <- c(sapply(x$protocols,function(u)u$treatment_variables))
    if (length(value$strategy) == 1 && value$strategy == "additive"){
        # FIXME: some formulas could be shared across protocols
        for (protocol in value$protocols){
            if (length(value[["exclude_other_treatments"]])>0){
                protocol_exclude_variables <- c(exclude_variables,setdiff(all_treatment_variables,x$protocols[[protocol]]$treatment_variables))
            } else{
                protocol_exclude_variables <- exclude_variables
            }
            x$models[[protocol]] = additive_formalizer(x = x,
                                                       protocol = protocol,
                                                       exclude_variables = protocol_exclude_variables,
                                                       include_variables = include_variables,
                                                       Markov = value$markov,
                                                       continuous_outcome = continuous_outcome)
            ## model(x) <- list(formalizer = "additive",treatment_variables = x$protocols[[value$protocol]]$treatment_variables)
            x$targets[[value$name]][["strategy"]] <- "additive"
        }
    }else{
        stop("Don't know about this strategy.")
    }
    x$targets[[value$name]][["protocols"]] <- unique(c(x$targets[[value$name]][["protocols"]],value$protocols))
    x$targets[[value$name]][["estimator"]] <- value$estimator
    x
}
######################################################################
### target.R ends here
