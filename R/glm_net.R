### glm_net.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 27 2024 (18:30) 
## Version: 
## Last-Updated: apr 23 2026 (08:26) 
##           By: Thomas Alexander Gerds
##     Update #: 15
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Fit glm_net models for use with rtmle
#'
#' @description Fits \code{glm_net} models from a formula and data set for use
#' with \code{\link{learn_glmnet}}.
#' @name glm_net
#'
#' @param formula A formula.
#' @param data The data on which to fit the model. 
#' @param lambda Tuning parameter for \code{glm_net}. If \code{NULL}, the
#'   parameter is chosen by \code{glmnet}.
#' @param cv Logical. If \code{TRUE}, use cross-validation.
#' @param alpha Elastic-net mixing parameter. See
#'   \code{\link[glmnet]{glmnet}}.
#' @param nfolds Number of folds for cross-validation. Default is 10.
#' @param type.measure Loss function used for cross-validation. Default is
#'   deviance.
#' @param family Passed to \code{\link[glmnet]{glmnet}}. Defaults to
#'   \code{"binomial"} for binary outcomes and \code{"gaussian"} otherwise.
#' @param ... Additional arguments passed to \code{\link[glmnet]{glmnet}}.
#' @return A fitted object of class \code{"glm_net"}.
#' @seealso \code{\link{learn_glmnet}}, \code{\link{run_rtmle}}
#' @examples
#' d <- data.frame(Y = rep(c(0, 1), 10),
#'                 A = rep(c(0, 1, 1, 0), 5),
#'                 L = seq(-1, 1, length.out = 20))
#' fit <- glm_net(Y ~ A + L, data = d, lambda = 0.01, cv = FALSE)
#' riskRegression::predictRisk(fit, newdata = d[1:3, ], times = 1, lambda = 0.01)
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
    if (is.character(formula)) formula <- formula(formula)
    tt <- all.vars(stats::update(formula,".~1"))
    sorted_x_train=bl_obj=terms=design = NULL
    y  <- data[[tt[1]]]
    if (missing(family)) {
        if (length(unique(y)) == 2){
            family <- "binomial"
        }else{
            family <- "gaussian"
        }
    }
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
