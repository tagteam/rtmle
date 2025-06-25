### print.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (08:31) 
## Version: 
## Last-Updated: Jun 25 2025 (12:55) 
##           By: Thomas Alexander Gerds
##     Update #: 57
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Printing an object for register data analysis with the targeted minimum loss estimator 
#'
#' This print function informs about the current state of the object and
#' about how to proceed.
#' @param x Object to be printed
#' @param ... Not used for now
#' @method print rtmle
#' @export
print.rtmle <- function(x,...){
    cat("\nTargeted minimum loss based analysis of longitudinal data on a discretized time scale",sep = "")
    cat("\ninitialized to ",length(x$times)-1," time intervals (starting at time zero).\n\n",sep = "")
    cat("See name of the subject id variable is '",x$names$id,"'.\nSee x$names for the names of the other variables.\n\n",sep = "")
    # data 
    has_data <- (length(x$data)>0|(has_long_data <- length(x$long_data)>0))
    if (!has_data)
        cat("\nTODO: The object contains no data yet.\n",
            "      - Add baseline data with the function 'add_baseline_data'\n",
            "      - Use the functions 'add_long_data' and 'long_to_wide' to add and discretize the dates of outcome, censored and competing events and to the longitudinal data (covariates and treatment).\n",
            "      - Alternatively, use the function 'add_wide_data' to add the readily discretized data.",
            sep = "")
    else{
        for (name in c("baseline_data","timevar_data","outcome_data")){
            if (length(x$data[[name]]) == 0 && length(x$long_data[[name]]) == 0) {
                cat("TODO: The object contains no ",
                    name,
                    ". Add them with the function ",
                    ifelse(name == "baseline_data","add_baseline_data", ifelse(has_long_data,"add_long_data","add_wide_data")),
                    ".\n",
                    sep = "")
            }
        }
    }
    # target trial protocols
    if (length(x$protocols)>0) {
        cat("The protocols are stored here: x$prototols. Add more protocols with the function 'protocol'.\n")
    } else{
        cat("TODO: The object contains no protocols. Add them with the function 'protocol'.\n",sep = "")        
    }
    if (length(x$targets)>0) {
        cat("The targets are stored here: x$targets\n")
        if (length(x$models) == 0){
            cat("TODO: Use the function 'model_formula' to initialize the formula for the nuisance parameter models.\n")
        }
    } else{
        cat("TODO: The object contains no targets yet. Add them with the function 'target'.\n",sep = "")        
    }
}


######################################################################
### print.rtmle.R ends here
