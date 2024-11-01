### update_Q.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:54) 
## Version: 
## Last-Updated: Nov  1 2024 (07:25) 
##           By: Thomas Alexander Gerds
##     Update #: 29
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
update_Q <- function(Y,
                     logitQ,
                     cum.g, 
                     uncensored_undeterministic,
                     intervention.match) {
    N <- length(Y)
    off <- as.vector(logitQ)
    if (length(cum.g) == 0) cum.g <- 1
    ## FIXME: are there better ways to remove those censored in current interval?
    subjects_with_weights <- !is.na(Y) & uncensored_undeterministic & as.vector(intervention.match)
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
        ## browser(skipCalls=TRUE)
        ## m <- ltmle.glm(f, data = data.temp[weights > 0, ], family = quasibinomial(),
        ## weights = as.vector(scale(weights[weights > 0], center = FALSE)))
        ## browser(skipCalls = TRUE)
        Qstar <- stats::predict(m, newdata = data.temp, type = "response")
    } else {
        warning("No TMLE update because no subject has positive weight")
        Qstar <- stats::plogis(logitQ)
        m <- "no Qstar fit because no subjects alive, uncensored, following intervention"
    }
    return(Qstar)
}



######################################################################
### update_Q.R ends here
