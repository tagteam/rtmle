### learn_glm.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 23 2024 (12:49) 
## Version: 
## Last-Updated: maj  7 2026 (09:45) 
##           By: Thomas Alexander Gerds
##     Update #: 154
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Learn nuisance-parameter models with glm
##'
##' Learns nuisance-parameter models for TMLE and predicts probabilities in
##' intervention-updated data using \code{\link[stats]{glm}}. The function first
##' attempts to use \code{\link[speedglm]{speedglm.wfit}} when available and
##' falls back to \code{\link[stats]{glm.fit}} if needed.
##'
##' @title Nuisance-parameter learner based on glm
##' @param character_formula Formula for the nuisance parameter, supplied as a
##'   character string.
##' @param data Data used for learning.
##' @param intervened_data Data used for prediction after intervention variables
##'   have been set according to a protocol.
##' @param learn_variables Optional vector of variable names to include in the
##'   learner. Use \code{"NONE"} to fit an intercept-only model that predicts
##'   the mean outcome.
##' @param maxit Number of iterations passed to \code{\link[stats]{glm.fit}}
##'   and \code{\link[speedglm]{speedglm.wfit}}.
#' @param save_fitted_objects Logical. If \code{TRUE}, store the 
#'   fitted object as element \code{fit}
##' @param reuse_fit Optional fitted object returned by a previous call. If
##'   supplied, skip fitting and only predict in \code{intervened_data}.
##' @param ... Additional arguments for the learning phase. Not currently used.
##' @return A list whose first element, \code{predicted_values}, is a vector of
##'   predicted probabilities. Element \code{fit} contains a coefficient summary,
##'   and element \code{object}, when available, contains the fitted model.
##' @seealso \code{\link{superlearn}}, \code{\link{learn_ranger}},
##'   \code{\link{learn_glmnet}}, \code{\link{learn_xgboost}}
##' @examples
##' d <- data.table::data.table(Y = rep(c(0, 1), 10),
##'                             A = rep(c(0, 1, 1, 0), 5),
##'                             L = seq(-1, 1, length.out = 20))
##' predicted <- learn_glm("Y ~ A + L", data = d, intervened_data = d)
##' head(predicted$predicted_values)
##' predicted$fit
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
learn_glm <- function(character_formula,
                      data,
                      intervened_data,
                      learn_variables = NULL,
                      maxit = 100,
                      save_fitted_objects = FALSE,
                      reuse_fit = NULL,
                      ...){
    # NOTE: speedglm could be handled by own learner 
    ## args <- as.list(match.call())
    ## if (!("speed"%in%names(args))) speed <- FALSE
    if (length(reuse_fit) == 0){    
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
                # learner_output
                fit = matrix(intercept,ncol = 1,nrow = 1,dimnames = list("(Intercept)","Estimate"))
                class(fit) <- c("constant_probability","list")
                predicted_values <- list(predicted_values = predicted_values,
                                         fit = fit,
                                         fit_summary = mean_Y)
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
        if (speed &&
            requireNamespace("speedglm", quietly = TRUE) &&
            !inherits(try(
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
        fit_summary <- coef(summary(fit))
    }else{
        fit <- reuse_fit$fit
        fit_summary = reuse_fit$fit_summary
    }
    predicted_values <- predict(
        fit,
        type = "response",
        newdata = intervened_data,
        se = FALSE
    )
    #parse_learner_output
    output <- list(predicted_values = predicted_values,
                   fit_summary = fit_summary)
    if (isTRUE(save_fitted_objects)){
        output$fit <- fit
    }
    output
}


######################################################################
### learn_glm.R ends here
