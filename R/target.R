### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: Sep 24 2024 (10:20) 
##           By: Thomas Alexander Gerds
##     Update #: 24
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
    if (length(value$strategy) == 1 && value$strategy == "additive"){
        # FIXME: some formulas could be shared across protocols
        for (protocol in value$protocols){
            x$models[[protocol]] = additive_formalizer(x = x,
                                                       protocol = protocol,
                                                       Markov = NULL)
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
