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
#' Adds long-format outcome, censoring, competing-event, treatment, and
#' time-varying covariate data.
##' @title Add long-format data
#' @param x An object of class \code{"rtmle"}, created by
#'   \code{\link{rtmle_init}}.
#' @param outcome_data A data frame with the subject identifier named as
#'   \code{x$names$id} and a column with outcome event times. Only subjects
#'   who experience an outcome event need rows.
#' @param censored_data A data frame with subject identifiers and
#'   end-of-follow-up times for subjects who are censored.
#' @param competing_data A data frame with subject identifiers and competing
#'   event times. Only subjects who experience a competing event need rows.
#' @param timevar_data A named list of data frames for treatment and
#'   time-varying covariates. Each data frame must contain the subject
#'   identifier and either a time column or the columns \code{start_date}
#'   and \code{end_date}. It may also contain a \code{value} column.
#'   Multiple rows per subject are allowed.
#' @details For discrete-time analyses, name the time column \code{date}.
#'   Supply calendar dates or numeric durations since the start of follow-up.
#'   If an outcome, censoring, or competing-event table contains multiple
#'   rows for one subject, the earliest row is stored.
#'
#'   For continuous-time analyses, use the time column configured in
#'   \code{\link{rtmle_init}} or \code{date}. Supply numeric durations
#'   since baseline and at most one row per subject in each outcome, censoring,
#'   or competing-event table. See \code{\link{add_continuous_data}} for
#'   continuous-time input and event definitions.
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
    id_name <- x$names$id
    continuous <- isTRUE(x$continuous_time)
    time_names <- if (continuous) unique(c(x$names$time, "date")) else "date"
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
    # values. Setting the following variables forces discretize_data() to
    # re-evaluate the date representation and the event histories.
    x$progress$substracted_start_followup_date <- FALSE
    x$progress$discretized_outcome <- FALSE
    x$progress$discretized_grid <- NULL
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
            if (!(id_name %in% names(current_data))){
                stop(paste0("Element ",
                            name,
                            " of argument 'timevar_data' does not have a variable called '",
                            id_name, "'."))
            }

            time_name <- intersect(time_names, names(current_data))[1L]
            has_interval <- all(c("start_date", "end_date") %in%
                                names(current_data))
            allowed_names <- if (has_interval) {
                c(id_name, "start_date", "end_date", "value")
            } else if (!is.na(time_name)) {
                c(id_name, time_name, "value")
            } else {
                character()
            }
            if (length(allowed_names) == 0L ||
                length(setdiff(names(current_data), allowed_names)) > 0L) {
                stop(paste0(
                    "Element ", name,
                    " of argument 'timevar_data' is inconsistent with the ",
                    "accepted formats: (id,", time_names[1L], "), (id,",
                    time_names[1L], ",value), ",
                    "(id,start_date,end_date), or ",
                    "(id,start_date,end_date,value)."
                ))
            }

            if (!continuous) {
                order_names <- c(id_name, if (has_interval) "start_date" else time_name)
                data.table::setorderv(current_data, order_names)
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
        time_name <- intersect(time_names, names(current_data))[1L]
        if (!id_name %in% names(current_data) || is.na(time_name)) {
            stop("Element ", name, " must contain variables called '", id_name,
                 "' and '", paste(time_names, collapse = "' or '"), "'.")
        }
        if (continuous) {
            if (anyDuplicated(current_data[[id_name]])) {
                stop("Duplicated values of variable '", id_name,
                     "' in dataset '", name, "' are not allowed")
            }
        } else {
            data.table::setorderv(current_data, c(id_name, time_name))
            current_data <- current_data[!duplicated(current_data[[id_name]])]
        }
        x$long_data[[name]] <- current_data
    }
    x
}


######################################################################
### add_long_data.R ends here
