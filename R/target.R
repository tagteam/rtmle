### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: dec  3 2025 (08:42) 
##           By: Thomas Alexander Gerds
##     Update #: 41
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
## ### Code:
##' Define a target parameter and specify formulas for the estimators of the nuisance parameters
##'
##' A target defines a target parameter, an estimator and models
##' for the nuisance parameters.
##' @param x An rtmle object as obtained with \code{rtmle_init}.
##' @param name the name of the target parameter.
##' @param protocols The names of the protocol(s) involved in the target parameters.
##' @param estimator Character specifying the estimator: either \code{'tmle'} or \code{'g-formula'}.
##' @param ... Not (yet) used
##' @export
target <- function(x,
                   name,
                   protocols,
                   estimator,
                   ...) {
    if (all(!(found <- protocols%in%names(x$protocols)))) {
        stop("The following protocols are not defined:\n",paste0(protocols[!found],collapse = ", "))
    }
    if (!missing(estimator)) warning("Argument estimator is obsolete and ignored. Specify new argument estimator of run_rtmle instead.")
    x$targets[[name]][["protocols"]] <- unique(c(x$targets[[name]][["protocols"]],protocols))
    x
}

######################################################################
### target.R ends here
