### add_long_data.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 25 2024 (11:24)
## Version:
## Last-Updated: Aug 12 2025 (19:25) 
##           By: Thomas Alexander Gerds
##     Update #: 51
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
## ### Code:
#' Adding long format data to a rtmle object
#'
#' This function adds a list of datasets with dates and values in long format (multiple lines per subject) to an existing rtmle object
##' @title Adding long format data
#' @param x object of class \code{rtmle}
#' @param outcome_data data.frame with two columns where the first is a subject identifier as initialized
#'        under \code{x$names$id} and the second column contains dates of outcome events. 
#'        Only subjects who exerience an outcome may have a row in the data.frame. The 
#'        name of the date variable must match \code{x$names$time} and it may contain calendar time dates
#'        or the pre-calculated time duration between time zero and the onset of the outcome. 
#' @param censored_data data.frame with two columns where the first is a subject identifier as initialized
#'        under \code{x$names$id} and the second column contains dates of end of followup.
#'        All subjects can have a row in the data.frame. But the entries about the end of followup for people who experience the outcome or 
#'        a competing risk are ignored.
#'        The name of the date variable must match \code{x$names$time} and it may contain calendar time dates
#'        or the pre-calculated time duration between time zero and the end of followup. 
#' @param competing_data data.frame with two columns where the first is a subject identifier as initialized
#'        under \code{x$names$id} and the second column contains dates of competing events. 
#'        Only subjects who exerience an competing event may have a row in the data.frame. The 
#'        name of the date variable must match \code{x$names$time} and it may contain calendar time dates
#'        or the pre-calculated time duration between time zero and the onset of the competing event. 
#' @param timevar_data named list of data.frames. 
#' @param ... Not used (not yet)
#' @seealso \link[rtmle]{add_baseline_data}, \link[rtmle]{add_wide_data}
#' @return The modified object.
#' @author Thomas A. Gerds <tag@@biostat.ku.dk>
#' @examples
#' set.seed(17)
#' tau <- 3
#' ld <- simulate_long_data(n = 91,number_visits = 20,
#'                          beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),
#'                          register_format = TRUE)
#' x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- add_long_data(x,
#'                    outcome_data=ld$outcome_data,
#'                    censored_data=ld$censored_data,
#'                    competing_data=ld$competing_data,
#'                    timevar_data=ld$timevar_data)
#' x <- add_baseline_data(x,data=ld$baseline_data)
#' x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
#' x <- protocol(x,name = "Always_A",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1"))))
#' x <- protocol(x,name = "Never_A",
#'                     intervention = data.frame("A" = factor("0",levels = c("0","1"))))
#' x <- prepare_data(x)
#' x <- target(x,name = "Outcome_risk",
#'                   estimator = "tmle",
#'                   protocols = c("Always_A","Never_A"))
#' x <- model_formula(x)
#' x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau)
#' summary(x)
#' @export
add_long_data <- function(x,
                          outcome_data,
                          censored_data,
                          competing_data,
                          timevar_data,
                          ...){
    nv = c(list(outcome_data,censored_data,competing_data))
    names(nv) = c("outcome_data","censored_data","competing_data")
    nv = nv[sapply(nv,NROW)>0]
    if (length(timevar_data)>0){
        if (!is.list(timevar_data) || is.data.frame(timevar_data) || is.null(names(timevar_data))) {
            stop("Argument timevar_data must be a named list of data.frames.")
        }
        for (name in names(timevar_data)){
            current_data <- copy(as.data.table(timevar_data[[name]]))
            if (!(all(c(x$names$id,"date") %in% names(current_data)))){
                warning(paste0("Element ",
                               name,
                               " of argument 'timevar_data' does not have variables called '",
                               x$names$id,
                               "' and 'date' and hence it is not added."))
            }else{
                if (length(names(current_data))>2) {
                    if (!(length(names(current_data)) == 3 & "value" %in%names(current_data)))
                        warning(paste0("Element ",
                                       name,
                                       " of argument 'timevar_data' has more than three variables or the 3rd variable is not called 'value' and hence it is not added."))
                }
                if (match(name,names(x$long_data$timevar_data),nomatch = 0)>0){
                    # replace existing element
                    x$long_data$timevar_data[[name]] <- current_data
                }
                else{
                    # append
                    tv <- list(current_data)
                    names(tv) <- name
                    x$long_data$timevar_data <- c(x$long_data$timevar_data,tv)
                }
            }
        }
    }
    for (name in names(nv)){
        if (!(x$names$id %in% names(nv[[name]])))
            warning(paste0("Element ",name," does not have a variable called ",x$names$id," and is not added."))
        else{
            x$long_data[[name]] <- copy(as.data.table(nv[[name]]))
        }
    }
    x
}


######################################################################
### add_long_data.R ends here
