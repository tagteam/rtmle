### rtmle_init.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (07:23) 
## Version: 
## Last-Updated: Oct  3 2024 (10:00) 
##           By: Thomas Alexander Gerds
##     Update #: 18
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Initialize an emulated trial analysis of register data for targeted minimum loss estimation  
##'
##' Create an empty object with class \code{"rtmle"} which waits for data
##' and other instructions such as   
##' @title Register Targeted Minimum Loss Estimation
##' @param intervals Number of discrete time intervals for updating information in a longitudinal setting.
##' Set this to 1 for settings where no updating of information should happen after time zero. 
##' @param name_id Name of the subject identifier variable
##' @param name_time Name of the time variable, e.g., \code{"date"} would be an appropriate name when the data are recorded as calendar dates.
##' @param name_outcome Name of the outcome variable(s). For example, the value \code{"cvddeath"} would mean
##' in a longitudinal setting with \code{intervals=3} that the variables in the data set are named \code{"cvddeath_1", "cvddeath_2", "cvddeath_3"}.
##' In the setting where \code{intervals=1} this should be the full name of the outcome variable.
##' @param name_competing Name of the competing risk variable(s).
##' @param name_censoring Name of the censoring variable(s).
##' @param censored_label Label of the values of the censoring variable(s) that indicated that
##' the data of the subject at this time interval are censored.
##' @return A list with class \code{"rtmle"} and the following elements:
##' \itemize{
##' \item
##' }
##' @seealso run_rtmle
##' @examples
##' # longitudinal setting
##' x <- rtmle_init(intervals=3,
##'                 name_id="pnr",
##'                 name_time="date",
##'                 name_outcome="cvddeath",
##'                 name_competing="death",
##'                 name_censoring="Censored",
##'                 censored_label="0")
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
rtmle_init <- function(intervals,
                       name_id,
                       name_time = "date",
                       name_outcome,
                       name_competing,
                       name_censoring = "Censored",
                       censored_label = "censored"){
    time_labels = paste0("time_",0:intervals)
    x = list(targets = NULL,
             estimate = NULL,
             names = list("id" = name_id,
                          "time" = name_time,
                          "outcome" = name_outcome,
                          "competing" = name_competing,
                          "censoring" = name_censoring,
                          "censored_label" = censored_label),
             times = 0:intervals)
    class(x) = "rtmle"
    x
}
######################################################################
### rtmle_init.R ends here
