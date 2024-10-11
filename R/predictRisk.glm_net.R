### predictRisk.glm_net.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 28 2024 (08:34) 
## Version: 
## Last-Updated: Oct  8 2024 (18:33) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## * predictRisk.glm_net
##' @export
##' @method predictRisk glm_net
predictRisk.glm_net <- function(object,newdata,times=NA,...) {
    args <- list(...)
    # check if user has supplied a lambda value
    slambda <- args$lambda
    # check if object has saved a selected lambda value
    if (length(slambda) == 0)
        slambda <- attr(object,"selected.lambda")
    if (length(slambda) == 0 || length(slambda)>1 || !(is.numeric(slambda))){
        stop("You must choose a single numeric lambda value for predictRisk ... ")
    }
    pos.lambda <- match(slambda,object$fit$lambda,nomatch = 0)
    if (pos.lambda == 0){
        stop("The fitted model was not fitted with the specified penalty parameter (lambda)")
    }
    lambda=cv=NULL
    # library(glmnet)
    # requireNamespace(c("prodlim","glmnet"))
    # predict.cv.glmnet <- utils::getFromNamespace("predict.cv.glmnet","glmnet")
    # predict.glmnet <- utils::getFromNamespace("predict.glmnet","glmnet")
    rhs <- stats::as.formula(stats::delete.response(object$terms))
    xnew <- stats::model.matrix(rhs,data=newdata)
    if (is.null(slambda) && object$cv){
        p <- stats::predict(object$fit,newx=xnew,type = "response", s="lambda.min")
    }
    else if (pos.lambda == 0 && !object$cv){
        if (length(object$lambda) == 1){
            p <- stats::predict(object$fit,newx=xnew,type = "response", s=object$lambda)
        }
        else {
            stop("Object fitted with multiple lambdas. You must pick one lambda for predictRisk!")
        }
    }
    else {
        p <- stats::predict(object$fit,newx=xnew,type = "response", s=slambda)
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    }
    p
}


######################################################################
### predictRisk.glm_net.R ends here
