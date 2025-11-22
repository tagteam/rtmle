### learn_ranger.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 28 2024 (09:26) 
## Version: 
## Last-Updated: nov 22 2025 (09:07) 
##           By: Thomas Alexander Gerds
##     Update #: 68
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Learning nuisance parameter models for TMLE and predicting
##' probabilities in intervened data based on \code{\link[ranger]{ranger}}
##'
##' Hyperparameters are set via ...
##' @title Nuisance parameter learner based on \code{\link[ranger]{ranger}}
##' @param character_formula Formula for nuisance parameter as a character
##' @param data Data for learning 
##' @param intervened_data Data for prediction 
##' @param ... Additional arguments for the learning phase passed to \code{\link[ranger]{ranger}}. These can include hyperparameters.
##' @return A vector of predicted probabilities which has the fit as an attribute.  
##' @seealso \code{link{superlearn}}, \code{link{learn_glm}}, \code{link{learn_glmnet}}
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
learn_ranger <- function(character_formula,data,intervened_data,...){
    # extract the data
    # FIXME: should not let NA's pass until here
    model_frame <- stats::model.frame(stats::formula(character_formula),
                                      data = data,
                                      drop.unused.levels = TRUE,
                                      na.action = na.omit)
    Y <- model_frame[[1]]
    if (length(unique(Y)) == 2) probability <- TRUE else probability <- FALSE
    if (probability) {
        if (!is.factor(Y))
            # FIXME: rtmle_predicted_outcome may only have two levels,
            #       and the values can be different from 0 and 1
            model_frame[[1]] <- factor(Y,levels = sort(unique(Y)))
    }
    if (!inherits(try(
             fit <- ranger::ranger(formula = stats::formula(character_formula),
                                   data = model_frame,
                                   probability = probability,
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
    no_missing <- !(apply(intervened_data,1,function(x)any(is.na(x))))
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
    data.table::setattr(predicted_values,"fit",NULL)
    return(predicted_values)
}


######################################################################
### learn_ranger.R ends here
