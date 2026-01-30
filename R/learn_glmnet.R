### learn_glmnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (16:42) 
## Version: 
## Last-Updated: jan 25 2026 (14:37) 
##           By: Thomas Alexander Gerds
##     Update #: 83
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Learning nuisance parameter models for TMLE and predicting
##' probabilities in intervened data based on \code{\link[glmnet]{glmnet}}
##'
##' This can be useful to avoid too many parameters in the model.
##' @title Nuisance parameter learner based on \code{\link[glmnet]{glmnet}}
##' @param character_formula Formula for nuisance parameter as a character
##' @param data Data for learning 
##' @param intervened_data Data for prediction
##' @param selector Character value deciding about how to select the penalty parameter lambda:
##' If \code{"undersmooth"} use the lambda value which results in the least amount of
##' penalty such that the model still fits. The other options are \code{"min"} and \code{"1se"} which are 
##' described in the documentation of \code{\link[glmnet]{cv.glmnet}}.
##'@param alpha The elasticnet mixing parameter. See the ?glmnet for more details.
#' @param nfolds Number of folds for cross-validation. Default is 10.
#' @param type.measure loss to use for cross-validation. Default is deviance.
##' @param ... Additional arguments for the learning phase passed to \code{\link[glmnet]{glmnet}}. 
##'        E.g., setting alpha affects the elastic net.
##' @return A vector of predicted probabilities which has the fit as an attribute.  
##' @seealso \code{link{superlearn}}, \code{link{learn_ranger}}, \code{link{learn_glm}}
#' @export
learn_glmnet <- function(character_formula,
                         data,
                         intervened_data,
                         selector = "undersmooth",
                         lambda = NULL,
                         alpha = 1,
                         nfolds = 10,
                         type.measure = "deviance",
                         ...){
    requireNamespace("glmnet")
    ## requireNamespace("riskRegression")
    RHS <- stats::formula(stats::delete.response(stats::terms(stats::formula(character_formula))))
    model_frame <- stats::model.frame(stats::formula(character_formula),
                                      data = data,
                                      drop.unused.levels = FALSE,
                                      na.action = na.omit)
    Y <- model_frame[[1]]
    ## FAM <- ifelse(length(unique(Y))>2,"gaussian","binomial")
    X <- model.matrix(RHS,data = data)
    if (selector == "undersmooth"){
        # Forcing cv = FALSE
        args <- c(list(cv = FALSE),list(x = X,y = Y,family = "gaussian",...))
        args <- args[unique(names(args))]
        ## fit <- do.call(glm_net,args)
        fit <- glmnet::glmnet(x=X,y=Y,lambda=lambda,alpha=alpha,...)
        selected.lambda <- fit$lambda[length(fit$lambda)]
        selected.beta <- fit$beta[,match(selected.lambda,fit$lambda,nomatch = NCOL(fit$beta)),drop = FALSE]
    }else{
        # forcing cv
        stopifnot(selector%in%c("min","1se"))
        args <- c(list(cv = TRUE),list(x = X,y = Y,family = "gaussian",...))
        args <- args[unique(names(args))]
        fit <- glmnet::cv.glmnet(x=X,y=Y,lambda=lambda,nfolds=nfolds,type.measure =type.measure,alpha=alpha,...)
        lambda <- fit$lambda
        selected.lambda <- fit[[paste0("lambda.",selector)]]
        selected.beta <- fit$glmnet.fit$beta[,fit$index[,"Lambda"][[selector]],drop = FALSE]
    }
    iX <- stats::model.matrix(RHS,data = intervened_data)
    predicted_values <- as.numeric(stats::predict(fit,newx=iX,type = "response", s=selected.lambda))
    data.table::setattr(predicted_values,"selected.lambda",selected.lambda)
    data.table::setattr(predicted_values,"fit",selected.beta)
    predicted_values
}


######################################################################
### learn_glmnet.R ends here
