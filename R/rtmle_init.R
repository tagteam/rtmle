### rtmle_init.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (07:23) 
## Version: 
## Last-Updated: apr 25 2026 (06:56) 
##           By: Thomas Alexander Gerds
##     Update #: 57
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
##' @param time_grid Vector of time points which define the discrete
##'     time intervals for updating information in a longitudinal
##'     setting. The first value must be 0 and the length of the
##'     vector must be at least 2. When \code{time_grid} is of length
##'     2, this is a setting where no updating of information happens
##'     after time zero.
##' @param name_id Name of the subject identifier variable
##' @param name_time Name of the time variable, e.g., \code{"date"} would be an appropriate name when the data are recorded as calendar dates.
##' @param name_outcome Name of the outcome variable(s). For example, the value \code{"cvddeath"} would mean
##' in a longitudinal setting with \code{time_grid} of length 4 that the variables in the data set are named \code{"cvddeath_1", "cvddeath_2", "cvddeath_3"}.
##' In the setting where \code{time_grid} has two values this should be the full name of the outcome variable.
##' @param name_competing Name of the competing risk variable(s). Can be \code{NULL} if the data do not contain competing risks for the outcome.
##' @param name_censoring Name of the censoring variable(s). Can be \code{NULL} if the subjects are all uncensored, i.e., the minimum
##' potential followup time is longer than the maximal time of \code{time_grid}.
##' @param censored_levels Character vector with the censoring levels.
##' @param censored_label A single character value. Label of the values of the censoring variable(s) that indicated that
##' the data of the subject at this time interval are censored. Must be an element of \code{censored_levels} too.
##' @param minority_threshold For binary outcomes, integer value which decides about whether
##' fitting nuisance parameter regression models is feasible. If the number of subjects in the minority group is
##' less than or equal to this threshold then regression modelling is skipped and the predicted value is simply the mean outcome.
##' @param weight_truncation Two values which decide about the minimum and the maximum of the weights w
##' used in inverse probability weighting (1/w).
##' @param prediction_range Two values used to force predicted outcome probabilities to the interval (0,1). This is required
##' to perform a logit transformation for the fit of the fluctuation model in the tmle update step. Default is \code{c(0.0001,0.9999)}.
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
##' x <- rtmle_init(time_grid=c(0,6,12,18),
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
rtmle_init <- function(time_grid,
                       name_id,
                       name_time = "date",
                       name_outcome,
                       name_competing = NULL,
                       name_censoring = NULL,
                       censored_levels = c("uncensored","censored"),
                       censored_label = "censored",
                       minority_threshold = 8,
                       weight_truncation = c(0,1),
                       prediction_range = c(0.0001,0.9999)){
    if(!(time_grid[1] == 0 & length(time_grid) > 1))
        stop("time_grid must be a vector of length > 1 starting with 0")
    if (length(name_censoring)>0){
        if (length(censored_label) != 1 ||
            length(censored_levels) != 2 ||
            !(censored_label %in% censored_levels) ||
            censored_label != censored_levels[[2]]){
            stop(paste0("Need exactly two values for censored_levels ",
                        "and one value for censored_label.\n",
                        "The censored_label must be the second value of censored_levels.")) 
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
             version = utils::packageVersion("rtmle"),
             time_grid_scale = time_grid,
             time_grid = 0:(length(time_grid)-1),
             intervention_nodes = 0:(length(time_grid)-2),
             tuning_parameters = list(minority_threshold = minority_threshold,
                                      weight_truncation = weight_truncation,
                                      prediction_range = prediction_range),
             diagnostics = NULL)
    class(x) = "rtmle"
    x
}
######################################################################
### rtmle_init.R ends here
