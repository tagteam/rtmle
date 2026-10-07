### add_long_data.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 25 2024 (11:24)
## Version:
## Last-Updated: okt  6 2026 (12:11) 
##           By: Thomas Alexander Gerds
##     Update #: 84
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
## ### Code:
#' Add long-format data to an rtmle object
#'
#' This function adds dates where outcome, censoring, treatment and confounders change
#' in long format to an rtmle object.
##' @title Add long-format data
#' @param x An object of class \code{"rtmle"}, created by
#'   \code{\link{rtmle_init}}.
#' @param outcome_data A data frame with two columns. The first column
#'   is the subject identifier and needs to be named as \code{x$names$id}.
#'   The second column is called \code{"date"} and contains outcome event times either as
#'   calendar dates or as numeric values. Date values are used by \code{discretize_data} to calculate
#'   durations from a subject specific start of followup date. Numeric values are typically  
#'   pre-computed durations since a start of followup date.
#'   Only subjects who experience an outcome event need rows. In case of multiple rows per
#'   subject only the first occurrence is stored.
#' @param censored_data Similar as \code{outcome_data} but with dates and duration indicating
#'   end-of-followup instead of outcome event. Only subjects who are lost to followup before the outcome/competing revent need rows.
#'   In case of multiple rows per
#'   subject only the first occurrence is stored.
#' @param competing_data Similar as \code{outcome_data} but with dates and duration indicating
#'   a competing event instead of an outcome event. Only subjects who experience a competing event need rows.
#' In case of multiple rows per
#'   subject only the first occurrence is stored.
#' @param timevar_data A named list of data frames for treatment and
#'   time-varying covariates. Each data frame must contain the subject
#'   identifier initialized as \code{x$names$id} and date information: 
#'   either a single variable called \code{"date"} or two columns called \code{"start_date"} and \code{"end_date"}. 
#'   The date information is either provided as calendar dates or as pre-computed durations since
#'   a start of followup date. In addition the data frame may 
#'   contain a \code{value} column. No entry, one row, or multiple rows per subject are allowed.
#' @param ... Named data-frame arguments for treatment or time-varying
#'   covariates. This is an alternative to passing those data frames through
#'   \code{timevar_data}.
#' @seealso \code{\link{rtmle_init}}, \code{\link{add_baseline_data}},
#'   \code{\link{add_wide_data}}, \code{\link{discretize_data}},
#'   \code{\link{prepare_rtmle_data}}
#' @return The modified \code{rtmle} object.
#' @author Thomas A. Gerds <tag@@biostat.ku.dk>
#' @examples
#' data(simulated_cohort)
#' ld <- register_format(simulated_cohort)
#' x <- rtmle_init(time_grid = seq(0,20,4),name_id = "id",
#'                 name_outcome = "stroke",name_competing = "death",
#'                 name_censoring = "dropout",censored_label = "censored")
#' x <- add_long_data(x,
#'                    outcome_data=ld$timevar_data$stroke[!duplicated(id)],
#'                    censored_data=ld$timevar_data$dropout,
#'                    competing_data=ld$timevar_data$death,
#'                    timevar_data=ld$timevar_data[c("bleeding","changeSBP","A","B")])
#' # alternative way to pass the time-varying covariates/treatments
#' x <- add_long_data(x,
#'                    outcome_data=ld$timevar_data$stroke[!duplicated(id)],
#'                    censored_data=ld$timevar_data$dropout,
#'                    competing_data=ld$timevar_data$death,
#'                    bleeding=ld$timevar_data[["bleeding"]],
#'                    changeSBP=ld$timevar_data[["changeSBP"]],
#'                    A=ld$timevar_data[["A"]],
#'                    B=ld$timevar_data[["B"]])
#' names(x$long_data)
#' x$long_data$outcome_data
#' x$long_data$competing_data
#' x$long_data$censored_data
#' x$long_data$timevar_data
#' @export
add_long_data <- function(x,
                          outcome_data,
                          censored_data,
                          competing_data,
                          timevar_data,
                          ...){
    if (is.null(x$progress)) x$progress <- list()
    if(missing(timevar_data)) timevar_data <- NULL
    if(missing(outcome_data)) outcome_data <- NULL
    if(missing(censored_data)) censored_data <- NULL
    if(missing(competing_data)) competing_data <- NULL
    if (length(timevar_data)>0){
        if (!is.list(timevar_data) || is.data.frame(timevar_data) || is.null(names(timevar_data))) {
            stop("Argument timevar_data must be a named list of data.frames.")
        }
    }
    # A new long-format data set may contain calendar dates even when a
    # previous call to discretize_data() readily converted dates to time-on-study
    # values. Setting the following variable, forces
    # discretize_data() to re-evaluate the date representation.
    x$progress$substracted_start_followup_date <- FALSE
    ## Add timevar_data supplied by dots:
    dots <- list(...)
    if (length(dots) == 0) {
        dots <- NULL
    } else {
        bad <- is.null(names(dots)) || any(names(dots) == "")
        if (bad) stop("All ... arguments must be named.")
    }
    timevar_data = c(timevar_data,dots)
    if(!is.null(timevar_data)){
        for (name in names(timevar_data)){
            if (!is.data.frame(timevar_data[[name]])) {
                stop("Element ", name,
                     " of argument 'timevar_data' must be a data.frame.")
            }
            current_data <- copy(as.data.table(timevar_data[[name]]))
            id_name <- x$names$id
            if (!(id_name %in% names(current_data))){
                stop(paste0("Element ",
                            name,
                            " of argument 'timevar_data' does not have a variable called '",
                            id_name, "'."))
            }

            has_date <- "date" %in% names(current_data)
            has_interval <- all(c("start_date", "end_date") %in%
                                names(current_data))
            allowed_names <- if (has_interval && !has_date) {
                c(id_name, "start_date", "end_date", "value")
            } else if (has_date && !has_interval) {
                c(id_name, "date", "value")
            } else {
                character()
            }
            if (length(allowed_names) == 0L ||
                length(setdiff(names(current_data), allowed_names)) > 0L) {
                stop(paste0(
                    "Element ", name,
                    " of argument 'timevar_data' is inconsistent with the " ,
                    "accepted formats: (id,date), (id,date,value), ",
                    "(id,start_date,end_date), or ",
                    "(id,start_date,end_date,value)."
                ))
            }

            order_names <- c(id_name,
                             if (has_interval) "start_date" else "date")
            data.table::setorderv(current_data, order_names)

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
    nv = c(list(outcome_data,censored_data,competing_data))
    names(nv) = c("outcome_data","censored_data","competing_data")
    # ignore empty data sets
    nv = nv[sapply(nv,NROW)>0]
    # check if id variable is stored
    for (name in names(nv)){
        if (!is.data.frame(nv[[name]])) {
            stop("Element ", name, " must be a data.frame.")
        }
        current_data <- copy(as.data.table(nv[[name]]))
        if (!(x$names$id %in% names(current_data)) ||
            !("date" %in% names(current_data))) {
            stop(paste0("Element ", name,
                        " must contain variables called '",
                        x$names$id, "' and 'date'."))
        }
        data.table::setorderv(current_data, c(x$names$id, "date"))
        current_data <- current_data[
            !duplicated(current_data[[x$names$id]])
        ]
        x$long_data[[name]] <- current_data
    }
    x
}


######################################################################
### add_long_data.R ends here
