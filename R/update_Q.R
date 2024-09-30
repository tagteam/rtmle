### update_Q.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:54) 
## Version: 
## Last-Updated: Sep 26 2024 (10:36) 
##           By: Thomas Alexander Gerds
##     Update #: 13
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
    subjects_with_weights <- uncensored_undeterministic & as.vector(intervention.match)
    weights <- numeric(N)
    weights[subjects_with_weights] <- 1/cum.g[subjects_with_weights]
    if (anyNA(weights)) stop("NA in weights")
    if (any(weights > 0)) {
        f <- as.formula("Y ~ -1 + S1 + offset(off)")
        data.temp <- data.frame(Y, S1 = rep(1,N),off)
        ## print(f)
        ## print(head(subjects_with_weights))
        ## print(tail(subjects_with_weights))
        ## print(head(as.vector(scale(weights[weights > 0], center = FALSE))))
        has_weight <- weights > 0
        weights <- as.vector(scale(weights[has_weight], center = FALSE))
        m <- glm(formula = f,
                 family = quasibinomial(),
                 data = data.frame(data.temp[has_weight, ],weights),
                 weights = weights,
                 control = glm.control(maxit = 100))
        ## browser(skipCalls=TRUE)
        ## m <- ltmle.glm(f, data = data.temp[weights > 0, ], family = quasibinomial(),
        ## weights = as.vector(scale(weights[weights > 0], center = FALSE)))
        ## browser(skipCalls = TRUE)
        Qstar <- predict(m, newdata = data.temp, type = "response")
    }
    else {
        Qstar <- plogis(logitQ)
        m <- "no Qstar fit because no subjects alive, uncensored, following intervention"
    }
    return(Qstar)
    ## indicator <- uncensored * intervention.match
    ## h.g.ratio <- indicator/cum.g 
    ## for (i in 1:num.regimes) {
    ## h.g.ratio[, i, ] <- h.g.ratio[, i, ] 
    ## weight.zero.index <- msm.weights[, i] == 0
    ## h.g.ratio[weight.zero.index, i, ] <- 0
    ## }
    ## return(list(Qstar = Qstar,
    ## h.g.ratio = h.g.ratio,
    ## off = off,
    ## fit = m))
}



######################################################################
### update_Q.R ends here
