### add_wide_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr  1 2025 (08:18) 
## Version: 
## Last-Updated: Apr  1 2025 (09:58) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Adding wide format data to a rtmle object
##'
#' This function adds a list of datasets in wide format (one line per subject) to an existing rtmle object
#' @title Adding wide format data
#' @param x object of class \code{rtmle} 
#' @param ... Not used (not yet)
#' @param value Named list of data.frames (or data.tables or tibbles). Two possible names 'outcome_data' and 'timevar_data'.
#' @return The modified object.
##' @seealso \link[rtmle]{add_baseline_data<-}, \link[rtmle]{add_long_data<-}
##' @examples
##' ## FIXME
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
"add_wide_data<-" <- function(x,...,value){
    nv <- tolower(names(value))
    if (!is.list(value) ||
        any(sapply(value,is.data.frame) == FALSE) ||
        any(is.na(nv)) ||
        any(nv == "") ||
        any(is.na(match(nv,c("outcome_data","timevar_data")))))
        stop(paste0("The 'value' must be a named list of data.frames (or data.tables or tibbles).\n",
                    "The names must match (one of):\n 'outcome_data','timevar_data'.\n"))
    if ("outcome_data"%in% nv){
        d <- data.table::copy(value[["outcome_data"]])
        if (!(x$names$id %in% names(d))){
            warning(paste0("Element 'outcome_data' does not have a variable called ",x$names$id," and is not added."))
        }else{
            if (!inherits(d,"data.frame")){
                stop("value[['outcome_data']] must be a data.frame (or data.table or tibble).")
            }else{
                ## check the registered names
                ## FIXME: could check if all Y_times are in the data
                ## FIXME: could check if all Y_0 is in the data which is not allowed
                for (Y in c("outcome","censoring","competing")){
                    if (length(x$names[[Y]])>0){
                        Y_variables <- grep(paste0("^",x$names[[Y]],"_[0-9]+"),names(d),value = TRUE)
                        if (length(Y_variables) == 0) {
                            stop(paste0("Cannot find the ",Y," variables: ",paste0(x$names[[Y]],"_",c("1","2","..."),collapse = ", ")))
                        }
                    }
                }
                x$data$outcome_data <- d
            }
        }
    }
    if ("timevar_data" %in% nv){
        varlist <- value[["timevar_data"]]
        # check if varlist is a list of data.frames
        if (!is.list(varlist) ||
            any(sapply(varlist,is.data.frame) == FALSE)){
            stop("value[['timevar_data']] must be a list of data.frames (or data.tables or tibbles).")
        }
        for (name in names(varlist)){
            d <- data.table::copy(varlist[[name]])
            if (!(x$names$id %in% names(d))){
                warning(paste0("Element 'timevar_data[['",name,"']] does not have a variable called ",x$names$id," and is not added."))
            } else{
                x$data$timevar_data[[name]] <- d
            }
        }
    }
    x
}


######################################################################
### add_wide_data.R ends here
