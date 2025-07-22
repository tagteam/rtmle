### make_exclusion_rule.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 22 2025 (09:06) 
## Version: 
## Last-Updated: Jul 22 2025 (09:18) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Helper function for model_formula
##'
##' This function creates a list of variables that should be excluded
##' from the formulas of a sequence of 
##' @title Helper function for model formula
##' @param timepoint Maximal time point of interest
##' @param name_outcome_variable Name of time-dependent outcome/treatment/censoring variable without suffix.
##' @param excluded_variables Vector of time-dependent variables (without suffix) to be excluded.
##' @param recursive Logical. If \code{TRUE} exclude all historical variables otherwise only the one from the same time point.
##' @return Named list 
##' @seealso model_formula
##' @examples
##' make_exclusion_rule(4,"A","B",recursive = TRUE)
##' make_exclusion_rule(4,"A",c("B","C"),recursive = TRUE)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
make_exclusion_rule <- function(timepoint,name_outcome_variable,excluded_variables,recursive = FALSE){
    timepoints <- 0:timepoint
    erule <- lapply(0:timepoint,function(tp){
        suffix <- if(recursive) 0:tp else tp
        unlist(lapply(excluded_variables,function(ev){paste0(ev,"_",suffix)}),use.names = FALSE)
    })
    names(erule) <- paste0(name_outcome_variable,"_",0:timepoint)
    erule
}


######################################################################
### make_exclusion_rule.R ends here
