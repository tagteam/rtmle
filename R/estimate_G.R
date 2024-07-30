### estimate_G.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (10:48) 
## Version: 
## Last-Updated: Jul 19 2024 (10:48) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

estimate_G <- function(x){
    n <- nrow(x$data)
    g <- cum.g <- prob.A.is.1 <- array(NaN, dim = c(n, length(x$Anodes)+length(x$Cnodes)))
    fit <- vector("list",length(x$Anodes)+length(x$Cnodes))
    names(fit) <- names(x$data)[nodes$AC]
    for (i in 1:(length(x$Anodes)+length(x$Cnodes))) {
        cur.node <- c(x$Anodes,x$Cnodes)[i]
        uncensored <- IsUncensored(x$uncensored,nodes$C,cur.node)
        deterministic.origdata <- IsDeterministic(x$data,
                                                  cur.node, x$deterministic.Q.function, nodes,
                                                  called.from.estimate.g = TRUE, x$survivalOutcome)$is.deterministic
        form <- x$gform[i]
        subs <- uncensored & !deterministic.origdata & !deterministic.g.origdata
        prob.A.is.1[, i, ] <- g.est$predicted.values
        fit[[i]] <- g.est$fit
        if (cur.node %in% nodes$A) {
            cur.abar <- AsMatrix(x$regimes[, nodes$A == cur.node, ])
            cur.abar.meanL <- cur.abar
        }
        else {
            cur.abar <- cur.abar.meanL <- matrix(1,nrow(x$data),1)
        }
        g[, i, ] <- CalcG(AsMatrix(prob.A.is.1[, i, ]),cur.abar,g.est$is.deterministic)
    }
    for (regime.index in 1:num.regimes) {
        cum.g.list <- CalcCumG(AsMatrix(g[, , regime.index]), x$gbounds)
        cum.g[, , regime.index] <- cum.g.list$bounded
    }
    return(list(cum.g = cum.g, fit = ReorderFits(fit), prob.A.is.1 = prob.A.is.1))
}

######################################################################
### estimate_G.R ends here
