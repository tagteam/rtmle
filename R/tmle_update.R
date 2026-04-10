### tmle_update.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:54) 
## Version: 
## Last-Updated: apr 10 2026 (15:22) 
##           By: Thomas Alexander Gerds
##     Update #: 62
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
                        intervention_probs, 
                        outcome_free_and_uncensored,
                        intervention_match,
                        k,
                        protocol) {
    N <- length(Y)
    if (length(intervention_probs) == 0) intervention_probs <- rep(1,N)
    ## FIXME: are there better ways to remove those censored in current interval?
    subjects_with_weights <- !is.na(Y) & outcome_free_and_uncensored & as.vector(intervention_match)
    weights <- numeric(N)
    weights[subjects_with_weights] <- 1/intervention_probs[subjects_with_weights]
    if (anyNA(weights)) stop("NA in weights")
    if (any(weights > 0)) {
        f <- stats::as.formula("Y ~ -1 + S1 + offset(offset)")
        data.temp <- data.frame(Y, S1 = rep(1,N),offset)
        has_weight <- weights > 0
        weights <- as.vector(scale(weights[has_weight], center = FALSE))
        m <- stats::glm(formula = f,
                        family = stats::quasibinomial(),
                        data = data.frame(data.temp[has_weight, ],weights),
                        weights = weights,
                        control = stats::glm.control(maxit = 100))
        Qstar <- stats::predict(m, newdata = data.temp, type = "response")
    } else {
        Qstar <- stats::plogis(offset)
        attr(Qstar,"diagnostics") <- data.table(Function = 'rtmle::tmle_update',
                                                Protocol = protocol,
                                                Step = k,
                                                Event = "No positive weights")
    }
    return(Qstar)
}



######################################################################
### tmle_update.R ends here
