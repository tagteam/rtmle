### print.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (08:31) 
## Version: 
## Last-Updated: Jul 29 2024 (14:51) 
##           By: Thomas Alexander Gerds
##     Update #: 13
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
#' @param ...
#' @method print rtmle
#' @export
print.rtmle <- function(x,...){
    cat("Targeted minimum loss analysis of register data.\n",sep = "")
    cat("  Initialized to ",length(x$times)," time intervals.\n",sep = "")
    cat("  The name of the subject id variable is '",x$name_id,"'.\n",sep = "")
    cat("  The outcome, competing risk and censoring variables are named '",x$name_outcome,"', '",x$name_competing,"', and '",x$name_censoring,"', respectively.\n",sep = "")
    for (name in c("baseline_data","timevar_data","treatment_data","outcome_data")){
        if (length(x$data[[name]])>0) {
        }else{
            cat("TODO: The object contains no ",name,". Add them with the function '",name,"<-'.\n",sep = "")
        }
    }
    if (length(x$targets)>0) {
    } else{
        cat("TODO: The object contains no protocols. Add them with the function 'protocol<-'.\n",sep = "")        
    }
    if (length(x$models)>0) {
    } else{
        cat("TODO: The object contains no models for the nuisance parameters. Add them with the function 'model<-'.\n",sep = "")        
    }
}


######################################################################
### print.rtmle.R ends here
