### learn_ranger.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 28 2024 (09:26) 
## Version: 
## Last-Updated: Oct 29 2024 (08:13) 
##           By: Thomas Alexander Gerds
##     Update #: 22
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
    if (!inherits(try(
             fit <- ranger::ranger(formula = stats::formula(character_formula),
                                   data = model_frame,
                                   probability = TRUE,
                                   ...)
            ,silent = TRUE),
             "try-error")){
    }else{
        stop(paste0("\nCould not fit model with ranger:\n",
                    "Formula:",character_formula))
    }
    predicted_values <- vector(mode = "numeric",length = NROW(intervened_data))
    this <- which(apply(intervened_data,1,function(x)any(is.na(x))))
    if (any(this)){
        predicted_values[this] <- as.numeric(NA)
        if (inherits(try(    
            predicted_values[!this] <- riskRegression::predictRisk(fit,newdata = intervened_data[!this])
        ),"try-error")) stop("Ranger prediction failed")
    }else{
        if (inherits(try(    
            predicted_values <- riskRegression::predictRisk(fit,newdata = intervened_data)
        ),"try-error")) stop("Ranger prediction failed")
    }
    ## print(predicted_values)
    data.table::setattr(predicted_values,"fit",NULL)
    return(predicted_values)
}


######################################################################
### learn_ranger.R ends here
