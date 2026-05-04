### learn_glmnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (16:42) 
## Version: 
## Last-Updated: feb 22 2026 (08:11) 
##           By: Thomas Alexander Gerds
##     Update #: 110
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Learn nuisance-parameter models with glmnet
##'
##' Learns nuisance-parameter models for TMLE and predicts probabilities in
##' intervention-updated data using \code{\link[glmnet]{glmnet}}.
##'
##' @title Nuisance-parameter learner based on glmnet
##' @param character_formula Formula for the nuisance parameter, supplied as a
##'   character string.
##' @param data Data used for learning.
##' @param intervened_data Data used for prediction after intervention variables
##'   have been set according to a protocol.
##' @param selector Character value controlling how the penalty parameter
##'   \code{lambda} is selected. If \code{"undersmooth"}, use the least
##'   penalized value that still fits the model. The other options,
##'   \code{"min"} and \code{"1se"}, are described in
##'   \code{\link[glmnet]{cv.glmnet}}.
##' @param lambda Penalty parameter passed to \code{\link[glmnet]{glmnet}} or
##'   \code{\link[glmnet]{cv.glmnet}}.
##' @param alpha Elastic-net mixing parameter. See
##'   \code{\link[glmnet]{glmnet}}.
##' @param family Passed to \code{\link[glmnet]{glmnet}} or
##'   \code{\link[glmnet]{cv.glmnet}}.
##' @param nfolds Number of folds for cross-validation. Default is 10.
##' @param type.measure Loss function used for cross-validation. Default is
##'   deviance.
##' @param ... Additional arguments passed to \code{\link[glmnet]{glmnet}}.
##' @return A list whose first element, \code{predicted_values}, is a vector of
##'   predicted probabilities. Element \code{fit} contains the selected
##'   coefficients, and element \code{object} contains the fitted model.
##' @seealso \code{\link{superlearn}}, \code{\link{learn_ranger}},
##'   \code{\link{learn_glm}}, \code{\link{learn_xgboost}}
##' @examples
##' d <- data.table::data.table(Y = rep(c(0, 1), 10),
##'                             A = rep(c(0, 1, 1, 0), 5),
##'                             L = seq(-1, 1, length.out = 20))
##' predicted <- learn_glmnet("Y ~ A + L", data = d, intervened_data = d,
##'                           lambda = 0.01)
##' head(predicted$predicted_values)
##' predicted$selected.lambda
#' @export
learn_glmnet <- function(character_formula,
                         data,
                         intervened_data,
                         selector = "undersmooth",
                         lambda = NULL,
                         alpha = 1,
                         family,
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
    Y <- as.numeric(model_frame[[1]])
    if (missing(family)){
        family <- ifelse(length(unique(Y))>2,"gaussian","binomial")
    }
    X <- model.matrix(RHS,data = data)
    args <- list(x = X,y = Y,lambda = lambda,alpha = alpha,family = family,...)
    args <- args[unique(names(args))]
    if (length(lambda) == 1 || selector == "undersmooth"){
        fit <- do.call(glmnet::glmnet,args)
        selected.lambda <- fit$lambda[length(fit$lambda)]
        selected.beta <- fit$beta[,match(selected.lambda,fit$lambda,nomatch = NCOL(fit$beta)),drop = FALSE]
    }else{
        # forcing cv
        stopifnot(selector%in%c("min","1se"))
        args <- c(list(type.measure = type.measure,nfolds = nfolds),args)
        args <- args[unique(names(args))]
        fit <- do.call(glmnet::cv.glmnet,args)
        ## lambda <- fit$lambda
        selected.lambda <- fit[[paste0("lambda.",selector)]]
        selected.beta <- fit$glmnet.fit$beta[,fit$index[,"Lambda"][[selector]],drop = FALSE]
    }
    iX <- stats::model.matrix(RHS,data = intervened_data)
    predicted_values <- as.numeric(stats::predict(fit,newx=iX,type = "response", s=selected.lambda))
    learner_output(predicted_values = predicted_values,
                   fit = selected.beta,
                   object = fit,
                   selected.lambda = selected.lambda)
}


######################################################################
### learn_glmnet.R ends here
