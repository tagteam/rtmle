### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Jul 26 2024 (11:54) 
##           By: Thomas Alexander Gerds
##     Update #: 10
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
    stopifnot(all(c("name","strategy","estimator","protocol")%in%names(value)))
    if (length(value$strategy) == 1 && value$strategy == "additive"){
        x$models = additive_formalizer(x = x,
                                       treatment_variables = x$protocols[[value$protocol]]$treatment_variables,
                                       Markov = NULL)
        ## model(x) <- list(formalizer = "additive",treatment_variables = x$protocols[[value$protocol]]$treatment_variables)
        x$targets[[value$name]][["strategy"]] <- "additive"
    }else{
        stop("Don't know about this strategy.")
    }
    x$targets[[value$name]][["protocol"]] <- value$protocol
    x$targets[[value$name]][["estimator"]] <- "tmle"
    x
}
######################################################################
### target.R ends here
