### add_wide_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr  1 2025 (08:18) 
## Version: 
## Last-Updated: Jul  1 2025 (11:38) 
##           By: Thomas Alexander Gerds
##     Update #: 26
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
#' @param outcome_data Data frame in a wide format. It should contain outcome variables,
#'             possibly also censoring and/or competing risk variables, and time-varying covariates.
#' @param timevar_data Data frame in a wide format. It should contain the time-varying covariates.
#' @param ... Not used (not yet)
#' @return The modified object.
##' @seealso \link[rtmle]{add_baseline_data}, \link[rtmle]{add_long_data}
##' @examples
#' # FIXME: this is a detour to produce wide format data
#' set.seed(112)
#' ld <- simulate_long_data(n = 11,number_visits = 20,
#'                          beta = list(A_on_Y = -.2,
#'                          A0_on_Y = -0.3,A0_on_A = 6),
#'                                register_format = TRUE)
#' x <- rtmle_init(intervals = 3, name_id = "id",
#'                 name_outcome = "Y", name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- add_long_data(x,
#'                    outcome_data=ld$outcome_data,
#'                    censored_data=ld$censored_data,
#'                    competing_data=ld$competing_data,
#'                    timevar_data=ld$timevar_data)
#' x <- add_baseline_data(x,data=ld$baseline_data)
#' x <- long_to_wide(x,intervals=seq(0,2000,30.45*6))
#' outcome_data <- x$data$outcome_data
#' timevar_data <- x$data$timevar_data
#' x <- rtmle_init(intervals = 3, name_id = "id",
#'                 name_outcome = "Y", name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- add_baseline_data(x,data=ld$baseline_data)
#' x <- add_wide_data(x,outcome_data=outcome_data,timevar_data=timevar_data)
#' x <- prepare_data(x)
#' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
add_wide_data <- function(x,outcome_data,timevar_data,...){
    if (!missing(outcome_data)){
        if (!inherits(outcome_data,"data.frame")){
            stop("outcome_data must be a data.frame (or data.table or tibble).")
        }
        if (!(x$names$id %in% names(outcome_data))){
            stop(paste0("Object 'outcome_data' does not have a variable called ",x$names$id,"."))
        }
        ## check the registered names
        ## FIXME: could check if all Y_times are in the data
        ## FIXME: could check if Y_0 is in the data which is not allowed
        for (Y in c("outcome","censoring","competing")){
            if (length(x$names[[Y]])>0){
                Y_variables <- grep(paste0("^",x$names[[Y]],"_[0-9]+$"),names(outcome_data),value = TRUE)
                if (length(Y_variables) == 0) {
                    stop(paste0("Cannot find the ",Y," variables: ",paste0(x$names[[Y]],"_",c("1","2","..."),collapse = ", ")))
                }
            }
        }
        x$data$outcome_data <- copy(outcome_data)
    }
    if (!missing(timevar_data)){
        if (is.data.frame(timevar_data)){
            d <- data.table::copy(timevar_data)
            if (!(x$names$id %in% names(d))){
                warning(paste0("Element 'timevar_data' does not have a variable called ",x$names$id," and is not added."))
            } else{
                x$data$timevar_data <- d
            }
        }else{
            x$data$timevar_data <- timevar_data            
        }
    }
    x
}


######################################################################
### add_wide_data.R ends here
