### learn_glmnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (16:42) 
## Version: 
## Last-Updated: Nov  4 2024 (06:58) 
##           By: Thomas Alexander Gerds
##     Update #: 36
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
learn_glmnet <- function(character_formula,data,intervened_data,selector = "undersmooth",...){
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
