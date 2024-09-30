### predictRisk.glm_net.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 28 2024 (08:34) 
## Version: 
## Last-Updated: Sep 28 2024 (08:35) 
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
## * predictRisk.glm_net
##' @rdname predictRisk
##' @method predictRisk glm_net
##' @export
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
    rhs <- as.formula(delete.response(object$terms))
    if (length(info <- object$surv_info) == 0){
        xnew <- model.matrix(rhs,data=newdata)
        if (is.null(slambda) && object$cv){
            p <- predict(object$fit,newx=xnew,type = "response", s="lambda.min")
        }
        else if (pos.lambda == 0 && !object$cv){
            if (length(object$lambda) == 1){
                p <- predict(object$fit,newx=xnew,type = "response", s=object$lambda)
            }
            else {
                stop("Object fitted with multiple lambdas. You must pick one lambda for predictRisk!")
            }
        }
        else {
            p <- predict(object$fit,newx=xnew,type = "response", s=slambda)
        }
    } else {
        # convert covariates to dummy variables
        newdata$dummy.time=rep(1,NROW(newdata))
        newdata$dummy.event=rep(1,NROW(newdata))
        dummy.formula=stats::update.formula(rhs,"prodlim::Hist(dummy.time,dummy.event)~.")
        EHF <- prodlim::EventHistory.frame(formula=dummy.formula,data=newdata,specials = NULL,unspecialsDesign=TRUE)
        newdata$dummy.time = NULL
        newdata$dummy.event = NULL
        # blank Cox object obtained with riskRegression:::coxModelFrame
        if (pos.lambda == 0 && object$cv){
            coxnet_pred <- c(exp(predict(object$fit,newx=EHF$design,type = "link", s="lambda.min")))
            lambda <- object$fit$lambda.min ## is needed for train_eXb
        }
        else if (pos.lambda == 0 && !object$cv){
            if (length(object$lambda) == 1){
                coxnet_pred <- c(exp(predict(object$fit,newx=EHF$design,type = "link", s=object$lambda)))
                lambda <- object$lambda
            }
            else {
                stop("Object fitted with multiple lambdas. You must pick a single value for lambda.")
            }
        } else {
            if (all((pos.lambda)>0)){
                coxnet_pred <- c(exp(predict(object$fit,newx=EHF$design,type = "link", s=slambda)))
                lambda <- slambda
            }
            else {
                stop(paste0("The fitted model was not fitted with the following penalty parameters (lambdas): ",
                            paste0(slambda[pos.lambda == 0],collapse = ", ")))
            }
        }
        train_eXb <- c(exp(predict(object$fit,newx=object$sorted_x_train,type = "link", s=lambda)))
        L0 <- riskRegression::baseHaz_cpp(starttimes = info$start,
                                          stoptimes = info$stop,
                                          status = info$status,
                                          eXb = train_eXb,
                                          strata = 1,
                                          nPatients = NROW(info$stop),
                                          nStrata = 1,
                                          emaxtimes = max(info$stop),
                                          predtimes = sort(unique(info$stop)),
                                          cause = 1,
                                          Efron = TRUE,
                                          reverse = FALSE)$cumhazard
        ## if (any(is.na(L0))) browser()
        coxnetSurv <- exp(-coxnet_pred%o%L0)
        where <- sindex(jump.times=unique(info$stop),eval.times=times)
        p <- cbind(0,1-coxnetSurv)[,1+where]
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    }
    p
}


######################################################################
### predictRisk.glm_net.R ends here
