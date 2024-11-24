### tmle_update.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:54) 
## Version: 
## Last-Updated: Nov 24 2024 (06:53) 
##           By: Thomas Alexander Gerds
##     Update #: 38
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
tmle_update <- function(Y,
                     offset,
                     cum.g, 
                     outcome_free_and_uncensored,
                     intervention.match) {
    N <- length(Y)
    if (length(cum.g) == 0) cum.g <- 1
    ## FIXME: are there better ways to remove those censored in current interval?
    subjects_with_weights <- !is.na(Y) & outcome_free_and_uncensored & as.vector(intervention.match)
    weights <- numeric(N)
    weights[subjects_with_weights] <- 1/cum.g[subjects_with_weights]
    if (anyNA(weights)) stop("NA in weights")
    if (any(weights > 0)) {
        f <- stats::as.formula("Y ~ -1 + S1 + offset(off)")
        data.temp <- data.frame(Y, S1 = rep(1,N),off)
        has_weight <- weights > 0
        weights <- as.vector(scale(weights[has_weight], center = FALSE))
        m <- stats::glm(formula = f,
                        family = stats::quasibinomial(),
                        data = data.frame(data.temp[has_weight, ],weights),
                        weights = weights,
                        control = stats::glm.control(maxit = 100))
        Qstar <- stats::predict(m, newdata = data.temp, type = "response")
    } else {
        warning("No TMLE update because no subject has positive weight")
        Qstar <- stats::plogis(offset)
        m <- "no Qstar fit because no subjects alive, uncensored, following intervention"
    }
    return(Qstar)
}



######################################################################
### tmle_update.R ends here
