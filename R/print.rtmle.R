### print.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (08:31) 
## Version: 
## Last-Updated: Apr  1 2025 (09:52) 
##           By: Thomas Alexander Gerds
##     Update #: 34
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
    cat("Targeted minimum loss analysis of register data.\n",sep = "")
    cat("  Initialized to ",length(x$times)-1," time intervals (starting at time zero).\n",sep = "")
    cat("  The name of the subject id variable is '",x$names$id,"'.\n",sep = "")
    cat("  The outcome, competing risk and censoring variables are named '",x$names$outcome,"', '",x$names$competing,"', and '",x$names$censoring,"', respectively.\n",sep = "")
    has_data <- length(x$data)>0
    if (!has_data)
        cat("TODO: The object contains no data yet.\nAdd the baseline data (if any) with the function 'add_baseline_data<-'\nAdd the outcome data and the longitudinal data (covariates and treatment) with one of the functions\n: 'add_long_data<-' or 'add_wide_data<-'.\n",
            sep = "")
    else{
        for (name in c("baseline_data","timevar_data","outcome_data")){
            if (length(x$data[[name]]) == 0 && length(x$long_data[[name]]) == 0) {
                cat("TODO: The object contains no ",name,". Add them with the function '",name,"<-'.\n",sep = "")
            }
        }
    }
    if (length(x$protocols)>0) {
    } else{
        cat("TODO: The object contains no protocols. Add them with the function 'protocol<-'.\n",sep = "")        
    }
    if (length(x$targets)>0) {
    } else{
        cat("TODO: The object contains no targets yet. Add them with the function 'target<-'.\n",sep = "")        
    }
}


######################################################################
### print.rtmle.R ends here
