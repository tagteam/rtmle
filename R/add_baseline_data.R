### add_baseline_data.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (08:54)
## Version:
## Last-Updated: Oct 14 2024 (07:13) 
##           By: Thomas Alexander Gerds
##     Update #: 9
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
### Code:
#' Adding baseline data to a rtmle object
#' 
#' This function adds a baseline dataset to an existing rtmle object
#' @param x object of class \code{rtmle} 
#' @param ... Not used (not yet)
#' @param value
#'@export
"add_baseline_data<-" <- function(x,...,value){
    if (!is.data.frame(value) || !(x$name$id %in% names(value)))
        stop(paste0("Value needs to be a data.frame (or data.table or tibble) with a variable called ",x$names$id))
    x$data$baseline_data <- value
    x
}



######################################################################
### add_baseline_data.R ends here
