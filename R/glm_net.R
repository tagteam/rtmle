### glm_net.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 27 2024 (18:30) 
## Version: 
## Last-Updated: Oct  2 2024 (16:02) 
##           By: Thomas Alexander Gerds
##     Update #: 10
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Fitting glm_net for use with rtmle
#'
#' @description Fit glm_net models via a formula and a data set for use with learn_glmnet
#' @name glm_net
#'
#' @param formula A formula.
#' @param data The data on which to fit the model. 
#' @param lambda The tuning parameters for glm_net. If set to NULL, then it the parameters are chosen for you.
#' @param cv Whether to use cross-validation or not. Default is TRUE.
#' @param alpha The elasticnet mixing parameter. See the ?glmnet for more details.
#' @param nfolds Number of folds for cross-validation. Default is 10.
#' @param type.measure loss to use for cross-validation. Default is deviance.
#' @param family passed to \code{glmnet}. Defaults for binary outcome to \code{"binomial"} and for survival to \code{"cox"}.
#' @param \dots Additional arguments that are passed on to the glmnet.
#' @export
glm_net <- function(formula,
                    data,
                    lambda=NULL,
                    cv=TRUE,
                    alpha = 1,
                    nfolds = 10,
                    type.measure = "deviance",
                    family,
                    ...){
    requireNamespace(c("glmnet","prodlim"))
    if (is.character(formula)) formula <- formula(formula)
    tt <- all.vars(stats::update(formula,".~1"))
    if (missing(family)) family <- "binomial"
    sorted_x_train=bl_obj=terms=design = NULL
    y  <- data[[tt[1]]]
    # factor levels should be ordered. e.g., first censored than uncensored
    if (is.factor(y)) {
        y <- as.numeric(y)-1
    }
    x <- stats::model.matrix(formula, data=data)
    if (!cv){
        fit <- glmnet::glmnet(x=x,
                              y=y,
                              lambda=lambda,
                              alpha=alpha,
                              family=family,
                              ...)
    }
    else {
        fit <- glmnet::cv.glmnet(x=x,
                                 y=y,
                                 lambda=lambda,
                                 nfolds=nfolds,
                                 type.measure =type.measure,
                                 alpha=alpha,
                                 family=family,
                                 ...)
        lambda <- fit$lambda
    }
    out = list(fit = fit,
               surv_info = bl_obj,
               design = design,
               call = match.call(),
               terms = terms(formula),
               sorted_x_train=sorted_x_train,
               cv=cv,
               lambda = lambda)
    class(out) = "glm_net"
    out
}


######################################################################
### glm_net.R ends here
