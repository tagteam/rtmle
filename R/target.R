### target.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 3 2024 (13:46)
## Version:
## Last-Updated: dec  4 2025 (12:04) 
##           By: Thomas Alexander Gerds
##     Update #: 43
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
## ### Code:
##' Define a target parameter
##'
##' A target records the regime(s) and estimator used for a target parameter.
##'
##' @param x An \code{rtmle} object as returned by \code{\link{rtmle_init}}.
##' @param name Name of the target parameter.
##' @param regimes Names of the regime(s) involved in the target parameter.
##' @param estimator Character string specifying the estimator: either
##'   \code{"tmle"} or \code{"g-formula"}.
##' @param ... Not used.
##' @return The modified \code{rtmle} object.
##' @seealso \code{\link{regime}}, \code{\link{model_formula}},
##'   \code{\link{run_rtmle}}, \code{\link{summary.rtmle}}
##' @examples
##' x <- rtmle_init(time_grid = 0:2, name_id = "id", name_outcome = "Y")
##' x <- regime(x, name = "Always_A",
##'               intervention = data.frame(time = x$intervention_nodes,
##'                                         A = factor("1", levels = c("0", "1"))))
##' x <- regime(x, name = "Never_A",
##'               intervention = data.frame(time = x$intervention_nodes,
##'                                         A = factor("0", levels = c("0", "1"))))
##' x <- target(x, name = "Outcome_risk",
##'             regimes = c("Always_A", "Never_A"),
##'             estimator = "tmle")
##' x$targets
##' x$estimator
##' @export
target <- function(x,
                   name,
                   regimes,
                   estimator,
                   ...) {
    if (all(!(found <- regimes%in%names(x$regimes)))) {
        stop("The following regimes are not defined:\n",paste0(regimes[!found],collapse = ", "))
    }
    x$targets[[name]][["regimes"]] <- unique(c(x$targets[[name]][["regimes"]],regimes))
    if (!missing(estimator)){
        x$estimator <- unique(c(x$estimator,estimator))
    }
    x
}

######################################################################
### target.R ends here
