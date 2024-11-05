### learn_ranger.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 28 2024 (09:26) 
## Version: 
## Last-Updated: Nov  5 2024 (13:41) 
##           By: Thomas Alexander Gerds
##     Update #: 41
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
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
            model_frame[[1]] <- factor(Y,levels = c(0,1))
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
                    "Formula:",character_formula))
    }
    predicted_values <- vector(mode = "numeric",length = NROW(intervened_data))
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
