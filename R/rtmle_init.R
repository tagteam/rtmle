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
##' Initialize an rtmle analysis
##'
##' Creates an empty object of class \code{"rtmle"} that stores analysis
##' settings and can then be populated with data, protocols, targets, and model
##' formulas.
##'
##' @title Register Targeted Minimum Loss Estimation
##' @param time_grid Vector of time points that define the discrete
##'     time intervals for updating information in a longitudinal
##'     setting. The first value must be 0 and the length of the
##'     vector must be at least 2. When \code{time_grid} is of length
##'     2, no information is updated
##'     after time zero.
##' @param name_id Name of the subject identifier variable.
##' @param name_time Name of the time variable. For example, \code{"date"} is
##'   appropriate when the data are recorded as calendar dates.
##' @param name_outcome Name of the outcome variable or variable prefix. For
##'   example, if \code{name_outcome = "cvddeath"} in a longitudinal setting
##'   with \code{time_grid} of length 4, the wide-format outcome variables are
##'   expected to be named \code{"cvddeath_1"}, \code{"cvddeath_2"}, and
##'   \code{"cvddeath_3"}. If \code{time_grid} has length 2, this should be
##'   the full name of the outcome variable.
##' @param name_competing Name of the competing-risk variable or variable
##'   prefix. Can be \code{NULL} if there are no competing risks for the
##'   outcome.
##' @param name_censoring Name of the censoring variable or variable prefix.
##'   Can be \code{NULL} if all subjects are uncensored; that is, if the
##'   minimum potential follow-up time is longer than the maximum value of
##'   \code{time_grid}.
##' @param censored_levels Character vector with the censoring levels.
##' @param censored_label A single character value identifying the censoring
##'   level that indicates censored observations. Must be an element of
##'   \code{censored_levels}.
##' @param minority_threshold For binary outcomes, an integer threshold used to
##'   decide whether fitting nuisance-parameter regression models is feasible.
##'   If the number of subjects in the minority group is less than or equal to
##'   this threshold, regression modeling is skipped and the predicted value is
##'   the mean outcome.
##' @param weight_truncation Two values defining the minimum and maximum
##'   inverse-probability weights.
##' @param prediction_range Two values used to constrain predicted outcome
##'   probabilities to the interval \code{(0, 1)}. This is required for the
##'   logit transformation used in the TMLE fluctuation model. The default is
##'   \code{c(0.0001, 0.9999)}.
##' @param time_grid_labels Labels used for the time-grid points on plot
##'   x-axes. Defaults to \code{time_grid}.
##' @return A list with class \code{"rtmle"} and the following elements:
##' \itemize{
##' \item targets
##' \item estimate
##' \item names
##' \item times
##' }
##' @seealso \code{\link{add_baseline_data}}, \code{\link{add_long_data}},
##'   \code{\link{add_wide_data}}, \code{\link{long_to_wide}},
##'   \code{\link{prepare_rtmle_data}}, \code{\link{protocol}},
##'   \code{\link{target}}, \code{\link{model_formula}},
##'   \code{\link{run_rtmle}}
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
                       prediction_range = c(0.0001,0.9999),
                       time_grid_labels = time_grid){
    if(!(time_grid[1] == 0 & length(time_grid) > 1))
        stop("time_grid must be a vector of length > 1 starting with 0")
    if (length(time_grid_labels) != length(time_grid)){
        stop("time_grid_labels must have the same length as time_grid")
    }
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
             time_grid_labels = as.character(time_grid_labels),
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
