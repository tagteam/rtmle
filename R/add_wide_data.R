### add_wide_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr  1 2025 (08:18) 
## Version: 
## Last-Updated: okt  6 2026 (14:15) 
##           By: Thomas Alexander Gerds
##     Update #: 48
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Add wide-format data to an rtmle object
##'
#' Adds one-row-per-subject data in wide format to an existing \code{rtmle}
#' object.
#'
#' @title Add wide-format data
#' @param x An object of class \code{"rtmle"}, typically created by
#'   \code{\link{rtmle_init}}.
#' @param outcome_data A wide-format data frame containing outcome variables and,
#'   when present, censoring and competing-risk variables.
#' @param timevar_data A wide-format data frame, or a named list of data frames,
#'   containing time-varying covariates.
#' @param ... Not used.
#' @return The modified \code{rtmle} object.
##' @seealso \code{\link{add_baseline_data}}, \code{\link{add_long_data}},
##'   \code{\link{discretize_data}}, \code{\link{prepare_rtmle_data}}
##' @examples
#' # Example data: Three subjects are observed at the time nodes 1, 2, and 3.
#' # Treatment and covariate values are available at the preceding nodes
#' # 0, 1, and 2.
#' outcome_long <- data.frame(
#'     id = rep(1:3, each = 3),
#'     time = rep(1:3, times = 3),
#'     Y = c(0, 1, 1,
#'           0, 0, 1,
#'           0, 0, 0)
#' )
#' treatment_long <- data.frame(
#'     id = rep(1:3, each = 3),
#'     time = rep(0:2, times = 3),
#'     A = c(0, 1, 1,
#'           1, 1, 0,
#'           0, 0, 1)
#' )
#' covariate_long <- data.frame(
#'     id = rep(1:3, each = 3),
#'     time = rep(0:2, times = 3),
#'     L = c(10, 12, 11,
#'           8, 9, 10,
#'           15, 14, 13)
#' )
#'
#' # Reshape each long-format table directly to the names expected by rtmle.
#' outcome_data <- reshape(outcome_long, idvar = "id", timevar = "time",
#'                         direction = "wide", sep = "_")
#' treatment_data <- reshape(treatment_long, idvar = "id", timevar = "time",
#'                           direction = "wide", sep = "_")
#' covariate_data <- reshape(covariate_long, idvar = "id", timevar = "time",
#'                           direction = "wide", sep = "_")
#'
#' x <- rtmle_init(time_grid = 0:3, name_id = "id", name_outcome = "Y")
#' x <- add_wide_data(
#'     x,
#'     outcome_data = outcome_data,
#'     timevar_data = list(A = treatment_data, L = covariate_data)
#' )
#' x$data$outcome_data
#' x$data$timevar_data
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
        for (Y in c("outcome","censoring","competing")){
            if (length(x$names[[Y]])>0){
                Y_variables <- grep(paste0("^",x$names[[Y]],"_[1-9]+$"),names(outcome_data),value = TRUE)
                zero_variable <- grep(paste0("^",x$names[[Y]],"_0$"),names(outcome_data),value = TRUE)
                if (length(zero_variable) > 0) {
                    if (length(unique(outcome_data[[zero_variable]]))>1){
                        stop(paste0("Found variation in the time zero variable ",zero_variable,".",
                                    "\nThe outcome, censoring, competing risk variables must all be zero at time zero."))
                    }
                }
                if (length(Y_variables) == 0) {
                    stop(paste0("Cannot find the ",Y," variables: ",paste0(x$names[[Y]],"_",c("1","2","..."),collapse = ", ")))
                }
            }
        }
        x$data$outcome_data <- data.table::copy(data.table::as.data.table(outcome_data))
        check_binary <- sapply(names(x$data$outcome_data),function(y){
            if ((y != x$names$id) &
                length(stats::na.omit(unique(x$data$outcome_data[[y]])))>2){
                stop(paste0("All outcome data must be binary. Variable '",y,"' has more than 2 values."))
            }
        })
    }
    if (!missing(timevar_data)){
        if ((!is.list(timevar_data) &&  !inherits(timevar_data,"data.frame")) || is.null(names(timevar_data))) {
            stop("Argument timevar_data must be a data.frame or a named list of data.frames.")
        }
        if (inherits(timevar_data,"data.frame")){
            d <- data.table::copy(data.table::as.data.table(timevar_data))
            if (!(x$names$id %in% names(d))){
                stop(paste0("Element 'timevar_data' does not have a variable called ",x$names$id,"."))
            } else{
                x$data$timevar_data <- d
            }
        }else{
            if (is.null(x$data$timevar_data)){
                x$data$timevar_data <- vector(length(timevar_data),mode = "list")
                names(x$data$timevar_data) <- names(timevar_data)
            }
            for (name in names(timevar_data)){
                current_data <- copy(as.data.table(timevar_data[[name]]))
                if (!(x$names$id %in% names(current_data))){
                    stop(paste0("Element 'timevar_data' does not have a variable called ",x$names$id,"."))
                } else{
                    x$data$timevar_data[[name]] <- current_data
                }
            }
        }
    }
    x
}


######################################################################
### add_wide_data.R ends here
