### learn_ranger.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 28 2024 (09:26) 
## Version: 
## Last-Updated: apr 29 2026 (07:32) 
##           By: Thomas Alexander Gerds
##     Update #: 81
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Learn nuisance-parameter models with ranger
##'
##' Learns nuisance-parameter models for TMLE and predicts probabilities in
##' intervention-updated data using \code{\link[ranger]{ranger}}.
##'
##' @title Nuisance-parameter learner based on ranger
##' @param character_formula Formula for the nuisance parameter, supplied as a
##'   character string.
##' @param data Data used for learning.
##' @param intervened_data Data used for prediction after intervention variables
##'   have been set according to a protocol.
##' @param ... Additional arguments passed to \code{\link[ranger]{ranger}},
##'   including hyperparameters.
##' @return A list whose first element, \code{predicted_values}, is a vector of
##'   predicted probabilities. Element \code{object} contains the fitted model.
##' @seealso \code{\link{superlearn}}, \code{\link{learn_glm}},
##'   \code{\link{learn_glmnet}}, \code{\link{learn_xgboost}}
##' @examples
##' d <- data.table::data.table(Y = factor(rep(c(0, 1), 10)),
##'                             A = rep(c(0, 1, 1, 0), 5),
##'                             L = seq(-1, 1, length.out = 20))
##' if (requireNamespace("ranger", quietly = TRUE)) {
##' predicted <- learn_ranger("Y ~ A + L", data = d, intervened_data = d,
##'                           num.trees = 10, min.node.size = 1)
##' head(predicted$predicted_values)
##' }
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
learn_ranger <- function(character_formula,data,intervened_data,...){
    if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("Package 'ranger' is required for learn_ranger().", call. = FALSE)
    }
    # extract the data
    model_frame <- stats::model.frame(stats::formula(character_formula),
                                      data = data,
                                      drop.unused.levels = TRUE,
                                      na.action = na.omit)
    Y <- as.numeric(model_frame[[1]])
    if (length(unique(Y)) == 2) probability <- TRUE else probability <- FALSE
    if (!inherits(
             try(
                 fit <- ranger::ranger(formula = stats::formula(character_formula),
                                       data = model_frame,
                                       probability = probability,
                                       write.forest = TRUE, importance = 'none', 
                                       ...)
                ,silent = TRUE),
             "try-error")){
    }else{
        stop(paste0("\nCould not fit model with ranger:\n",
                    "Formula: ",character_formula))
    }
    predicted_values <- vector(mode = "numeric",length = NROW(intervened_data))
    ivars <- all.vars(formula(character_formula))[-1]
    # remove all other variables to avoid false positive missing values
    intervened_data <- intervened_data[,ivars,with = FALSE]
    # then check for missing values
    no_missing <- !rowSums(is.na(intervened_data))
    predicted_values[!no_missing] <- as.numeric(NA)
    if (probability) {
        if (inherits(try(
            predicted_values[no_missing] <- riskRegression::predictRisk(fit,newdata = intervened_data[no_missing])
        ),"try-error")) {
            stop("Ranger prediction failed")
        }
    }else{
        if (inherits(try(
            predicted_values[no_missing] <- stats::predict(fit,data = intervened_data[no_missing],importance = "none")$predictions
        ),"try-error")) {
            stop("Ranger prediction failed")
        }
    }
    learner_output(predicted_values = predicted_values,
                   object = fit)
}


######################################################################
### learn_ranger.R ends here
