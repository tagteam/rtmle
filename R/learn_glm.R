### rtmle_glm.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (12:49) 
## Version: 
## Last-Updated: Nov  1 2024 (09:19) 
##           By: Thomas Alexander Gerds
##     Update #: 54
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
learn_glm <- function(character_formula,data,intervened_data,learn_variables = NULL,...){
    # FIXME: maxit should be controllable
    speed = TRUE
    # extract the data
    # FIXME: should not let NA's pass until here
    if (length(learn_variables)>0){
        allvars <- all.vars(stats::formula(character_formula))
        outcome_variable <- allvars[1]
        keepvars <- intersect(learn_variables,allvars[-1])
        if (length(keepvars)>0)
            character_formula <- paste0(outcome_variable,"~",paste0(keepvars,collapse = "+"))
        else
            character_formula <- paste0(outcome_variable,"~1")
    }
    model_frame <- stats::model.frame(stats::formula(character_formula),
                                      data = data,
                                      drop.unused.levels = TRUE,
                                      na.action = na.pass)
    Y <- model_frame[[1]]
    Y_label <- names(model_frame)[[1]]
    tf <- stats::terms(model_frame)
    X <- stats::model.matrix(object = tf, data = model_frame)
    if (any(is.na(X))) stop("Missing values in X")
    speed = TRUE
    if (speed && !inherits(try(
                      fit <- speedglm::speedglm.wfit(y = Y[!is.na(Y)],
                                                     X = X[!is.na(Y),],
                                                     intercept = attributes(tf)$intercept,
                                                     offset = stats::model.offset(model_frame),
                                                     family = quasibinomial(),
                                                     maxit = 100)
                     ,silent = TRUE),
                      "try-error")){
        class(fit) <- c("speedglm", "speedlm")
    }else{
        if (!inherits(try(
                 fit <- stats::glm.fit(y = Y[!is.na(Y)],
                                       x = X[!is.na(Y),],
                                       offset = stats::model.offset(model_frame),
                                       intercept = attributes(tf)$intercept,
                                       family = quasibinomial(),
                                       control = glm.control(maxit = 100))
                ,silent = TRUE),
                 "try-error")){
            class(fit) <- c(class(fit),"glm")
        }else{
            stop(paste0("\nCould not fit model with glm:\n",
                        "Outcome: ",
                        names(model_frame)[[1]],
                        "\nRight hand side: ",
                        paste0(colnames(X),collapse = "+ ")))
        }
        class(fit) <- c(class(fit),"glm")
    }
    fit$terms <- tf
    predicted_values <- predict(fit, type = "response", newdata = intervened_data, se = FALSE)
    ## print(predicted_values)
    data.table::setattr(predicted_values,"fit",coef(summary(fit)))
    return(predicted_values)
}


######################################################################
### rtmle_glm.R ends here
