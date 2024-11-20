### learn_glmnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (16:42) 
## Version: 
## Last-Updated: Nov 20 2024 (11:49) 
##           By: Thomas Alexander Gerds
##     Update #: 43
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Learning nuisance parameter models for TMLE and predicting
##' probabilities in intervened data based on \code{\link{glmnet}}
##'
##' This can be useful to avoid too many parameters in the model.
##' @title Nuisance parameter learner based on \code{\link{glmnet}}
##' @param character_formula Formula for nuisance parameter as a character
##' @param data Data for learning 
##' @param intervened_data Data for prediction
##' @param selector Character value deciding about how to select the penalty parameter lambda:
##' If \code{"undersmooth"} use the lambda value which results in the least amount of
##' penalty such that the model still fits. The other options are \code{"lambda.min"} and \code{"lambda.1se"} which are 
##' described in the documentation of \code{\link{glmnet}}.
##' @param ... Additional arguments for the learning phase passed to \code{\link{glmnet}}. 
##'        E.g., setting alpha affects the elastic net.
##' @return A vector of predicted probabilities which has the fit as an attribute.  
##' @seealso \code{link{superlearn}}, \code{link{learn_ranger}}, \code{link{learn_glm}}
#' @export
learn_glmnet <- function(character_formula,
                         data,
                         intervened_data,
                         selector = "undersmooth",
                         ...){
    requireNamespace("glmnet")
    ## requireNamespace("riskRegression")
    ## FAM <- ifelse(length(unique(Y))>2,"gaussian","binomial")
    model_frame <- stats::model.frame(stats::formula(character_formula),
                                      data = data,
                                      drop.unused.levels = FALSE,
                                      na.action = na.omit)
    if (selector == "undersmooth"){
        # Forcing cv = FALSE
        args <- c(list(cv = FALSE),list(formula = character_formula,
                                        data = model_frame,
                                        family = "gaussian",
                                        ...))
        args <- args[unique(names(args))]
        fit <- do.call(glm_net,args)
        selected.lambda <- fit$fit$lambda[length(fit$fit$lambda)]
    }else{
        # forcing cv
        stopifnot(selector%in%c("lambda.min","lambda.1se"))
        args <- c(list(cv = TRUE),list(formula = character_formula,data = data,family = "gaussian",...))
        args <- args[unique(names(args))]
        fit <- do.call(glm_net,args)
        selected.lambda <- fit$fit[[selector]]
    }
    imodel_frame <- stats::model.frame(stats::formula(stats::delete.response(stats::terms(stats::formula(character_formula)))),
                                       data = intervened_data,
                                       drop.unused.levels = FALSE,
                                       na.action = na.fail)
    predicted_values <- predictRisk(fit,
                                    type = "response",
                                    newdata = imodel_frame,
                                    lambda = selected.lambda)

    data.table::setattr(predicted_values,"selected.lambda",selected.lambda)
    predicted_values
}


######################################################################
### learn_glmnet.R ends here
