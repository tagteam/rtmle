### rtmle_init.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (07:23) 
## Version: 
## Last-Updated: Jun 16 2025 (14:17) 
##           By: Thomas Alexander Gerds
##     Update #: 38
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
##' @param name_competing Name of the competing risk variable(s). Can be \code{NULL} if the data do not contain competing risks for the outcome.
##' @param name_censoring Name of the censoring variable(s). Can be \code{NULL} if the subjects are all uncensored, i.e., the minimum
##' potential followup time is longer than the time of the last interval. 
##' @param censored_levels Character vector with the censoring levels.
##' @param censored_label A single character value. Label of the values of the censoring variable(s) that indicated that
##' the data of the subject at this time interval are censored. Must be an element of \code{censored_levels} too.
##' @return A list with class \code{"rtmle"} and the following elements:
##' \itemize{
##' \item targets
##' \item estimate
##' \item names
##' \item times
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
##'                 censored_levels = c("1","0"),
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
                       censored_levels = c("uncensored","censored"),
                       censored_label = "censored"){
    time_labels = paste0("time_",0:intervals)
    if (length(name_censoring)>0){
        if (length(censored_label) != 1 ||
            length(censored_levels) != 2 ||
            !(censored_label %in% censored_levels)){
            stop(paste0("Expect two values for censored_levels ",
                        "and one value for censored_label.\n",
                        "The censored_label must be an element of censored_levels.")) 
        }
        uncensored_label = setdiff(censored_levels,censored_label)
    }else{
        uncensored_label <- NULL
    }
    x = list(targets = NULL,
             estimate = NULL,
             data = NULL,
             names = list("id" = name_id,
                          "time" = name_time,
                          "outcome" = name_outcome,
                          "competing" = name_competing,
                          "censoring" = name_censoring,
                          ## "treatment_options" = treatment_options,
                          "censored_levels" = censored_levels,
                          "censored_label" = censored_label,
                          "uncensored_label" = uncensored_label),
             times = 0:intervals)
    class(x) = "rtmle"
    x
}
######################################################################
### rtmle_init.R ends here
