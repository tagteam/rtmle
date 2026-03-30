### learn_glm.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (12:49) 
## Version: 
## Last-Updated: mar 25 2026 (11:55) 
##           By: Thomas Alexander Gerds
##     Update #: 142
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Learning nuisance parameter models for TMLE and predicting
##' probabilities in intervened data based on \code{\link{glm}}
##'
##' The function attempts to run \code{speedglm} and uses \code{glm} only as a fallback option.
##' Optional argument learn_variables can be used to restrict the learning to a subset of variables.
##' This can be useful to avoid too many parameters in the model.
##' @title Nuisance parameter learner based on \code{\link{glm}}
##' @param character_formula Formula for nuisance parameter as a character
##' @param data Data for learning 
##' @param intervened_data Data for prediction 
##' @param learn_variables Vector of variable names. Only include these in the learning.
##' Can be \code{"NONE"} which means to not use any variables and simply predict the mean outcome.
##' @param maxit Number of iterations passed to \code{link[stats]{glm.fit}} and \code{link[speedglm]{speedglm.wfit}} 
##' @param ... Additional arguments for the learning phase. Not used at the moment.
##' @return A vector of predicted probabilities which has the fit as an attribute.  
##' @seealso \code{link{superlearn}}, \code{link{learn_ranger}}, \code{link{learn_glmnet}}
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
learn_glm <- function(character_formula,
                      data,
                      intervened_data,
                      learn_variables = NULL,
                      maxit = 100,
                      ...){
    # FIXME: speedglm should be handled by own learner 
    ## args <- as.list(match.call())
    ## if (!("speed"%in%names(args))) speed <- FALSE
    speed <- TRUE
    # extract the data
    if (length(learn_variables)>0){
        allvars <- all.vars(stats::formula(character_formula))
        outcome_variable <- allvars[1]
        ## keepvars <- intersect(learn_variables,allvars[-1])
        keepvars <- grep(paste0("(^",learn_variables,"($|_[0-9][0-9]?$))",collapse = "|"),
                         allvars[-1],
                         value = TRUE)
        if (length(keepvars)>0){
            character_formula <- paste0(outcome_variable,"~",paste0(keepvars,collapse = "+"))
        } else{
            # No covariates
            Y <- data[[outcome_variable]]
            mean_Y <- mean(Y,na.rm = TRUE)
            predicted_values <- rep(mean_Y,NROW(intervened_data))
            intercept <- log(mean_Y/(1-mean_Y))
            data.table::setattr(predicted_values,"fit",matrix(intercept,ncol = 1,nrow = 1,dimnames = list("(Intercept)","Estimate")))
            return(predicted_values)
        }
    }
    model_frame <- stats::model.frame(stats::formula(character_formula),
                                      data = data,
                                      drop.unused.levels = TRUE,
                                      na.action = na.pass)
    Y <- as.numeric(model_frame[[1]])
    Y_label <- names(model_frame)[[1]]
    tf <- stats::terms(model_frame)
    X <- stats::model.matrix(object = tf, data = model_frame)
    ## FIXME: need to check missing values much earlier
    if (any(is.na(X))) {
        stop("Missing values in X")
    }
    if (speed && !inherits(try(
                      fit <- speedglm::speedglm.wfit(y = Y[!is.na(Y)],
                                                     X = X[!is.na(Y),],
                                                     intercept = attributes(tf)$intercept,
                                                     offset = stats::model.offset(model_frame),
                                                     family = quasibinomial(),
                                                     maxit = maxit),
                      silent = TRUE),
                      "try-error")){
        class(fit) <- c("speedglm", "speedlm")
    }else{
        if (!inherits(try(
                 fit <- stats::glm.fit(y = Y[!is.na(Y)],
                                       x = X[!is.na(Y),],
                                       offset = stats::model.offset(model_frame),
                                       intercept = attributes(tf)$intercept,
                                       family = quasibinomial(),
                                       control = glm.control(maxit = maxit))
                ,silent = TRUE),
                 "try-error")){
            class(fit) <- c(class(fit),c("glm","lm"))
        }else{
            stop(paste0("\nCould not fit model with glm:\n",
                        "Outcome: ",
                        names(model_frame)[[1]],
                        "\nRight hand side: ",
                        paste0(colnames(X),collapse = "+ ")))
        }
    }
    fit$terms <- tf
    # FIXME: need to collect these warnings instead of preventing them
    suppressWarnings(predicted_values <- predict(fit,type = "response",newdata = intervened_data,se = FALSE))
    data.table::setattr(predicted_values,"fit",coef(summary(fit)))
    return(predicted_values)
}


######################################################################
### learn_glm.R ends here
