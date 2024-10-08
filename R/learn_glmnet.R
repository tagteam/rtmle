### learn_glmnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (16:42) 
## Version: 
## Last-Updated: Sep 28 2024 (08:28) 
##           By: Thomas Alexander Gerds
##     Update #: 16
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
learn_glmnet <- function(formula,data,selector = "undersmooth",...){
    requireNamespace("glmnet")
    requireNamespace("riskRegression")
    ## FAM <- ifelse(length(unique(Y))>2,"gaussian","binomial")
    if (selector == "undersmooth"){
        # Forcing cv = FALSE
        args <- c(list(cv = FALSE),list(formula = formula,data = data,family = "gaussian",...))
        args <- args[unique(names(args))]
        fit <- do.call(glm_net,args)
        selected.lambda <- fit$fit$lambda[length(fit$fit$lambda)]
    }else{
        # forcing cv
        stopifnot(selector%in%c("lambda.min","lambda.1se"))
        args <- c(list(cv = TRUE),list(formula = formula,data = data,family = "gaussian",...))
        args <- args[unique(names(args))]
        fit <- do.call(glm_net,args)
        selected.lambda <- fit$fit[[selector]]
    }
    attr(fit,"selected.lambda") <- selected.lambda
    fit
}


######################################################################
### learn_glmnet.R ends here
