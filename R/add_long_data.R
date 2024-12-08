### add_long_data.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 25 2024 (11:24)
## Version:
## Last-Updated: Nov 25 2024 (08:04) 
##           By: Thomas Alexander Gerds
##     Update #: 16
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
## ### Code:
#' Adding register data to a rtmle object
#'
#' This function adds a list of datasets with dates and values in long format (multiple lines per subject) to an existing rtmle object
#' @param x object of class \code{rtmle} 
#' @param ... Not used (not yet)
#' @param value Named list of data.frames or data.tables or tibbles. Each data.frame
#' has to contain a subject identifier variable which must have the same name as defined
#' by \code{x$names$id}. However, not all subjects must have a row in the data.frame.
#' Each data.frame should also have a time or date variable which must have the same name
#' as defined by \code{x$names$time}.
#' @export
"add_long_data<-" <- function(x,...,value){
    nv <- names(value)
    if (!is.list(value) || any(sapply(value,is.data.frame) == FALSE) || any(is.na(nv)) || any(nv == ""))
        stop("value must be a named list of data.frames (or data.tables or tibbles)")
    for (name in nv){
        d <- value[[name]]
        if (!(x$names$id %in% names(d)))
            warning(paste0("Element ",name," lacks the identifier variable ",x$names$id," and is not added."))
        else{
            x$register_data <- copy(d)
        }
    }
    x
}


######################################################################
### add_long_data.R ends here
