### discretize_exposure.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 24 2026 (11:04) 
## Version: 
## Last-Updated: apr 10 2026 (15:46) 
##           By: Thomas Alexander Gerds
##     Update #: 135
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# data contains exposure episodes
# grid contains intervals on a discrete time scale
# the function discretizes how many days a subject
# was exposed within the interval and defines exposure when the
# number of days exceeds threshold. E.g., when the length of the
# interval is 6 months (180 days) and the threshold is set to 91 then exposure
# is 1 if the treatment duration exceeds 50% of the length interval
# and 0 otherwise
#' Discretize long-format time-dependent data onto a discrete time grid
#'
#' This function maps long-format event, measurement, or exposure histories
#' onto a discrete time grid and returns a wide-format representation with
#' one column per time interval.
#'
#' It is the core engine used by \code{long_to_wide()} and can also be used directly.
#'
#' @param method Character string specifying how the long-format data should
#'   be mapped to the discrete grid. Built-in methods include:
#'   \describe{
#'     \item{"measurement"}{Aggregate measurements within intervals.}
#'     \item{"locf"}{Last observation carried forward.}
#'     \item{"event"}{Event indicator carried forward after occurrence.}
#'     \item{"event_interval"}{Indicator for events occurring within an interval.}
#'     \item{"exposure_time"}{Total exposure time within interval.}
#'     \item{"exposure_percent"}{Fraction of interval exposed.}
#'     \item{"any_exposure"}{Indicator for any exposure above threshold.}
#'     \item{"has_exposure"}{Indicator for exposure exceeding threshold proportion.}
#'   }
#'
#' @param data A \code{data.table} containing long-format observations.
#'   Required columns depend on \code{method}:
#'   \describe{
#'     \item{event / locf}{\code{id}, \code{date}}
#'     \item{measurement}{\code{id}, \code{date}, \code{value}}
#'     \item{exposure_*}{\code{id}, \code{start_date}, \code{end_date}}
#'   }
#'
#' @param grid A \code{data.table} defining the discrete time grid with
#'   one row per subject and interval. Must contain:
#'   \describe{
#'     \item{\code{id}}{Subject identifier}
#'     \item{\code{interval}}{Discrete time index}
#'     \item{\code{start_interval}, \code{end_interval}}{Interval boundaries}
#'   }
#'
#' @param name Character string used as prefix for the resulting wide-format
#'   variables.
#'
#' @param id Character string naming the subject identifier column.
#'
#' @param threshold Numeric threshold used for methods such as
#'   \code{"any_exposure"} and \code{"has_exposure"}.
#'
#' @param lookback_window Numeric or \code{Inf}. For methods using rolling joins
#'   (e.g., \code{"locf"} and \code{"event"}), defines how far back in time
#'   values are carried forward.
#'
#' @param values Length-2 vector specifying values used for binary indicators,
#'   typically \code{c(1, 0)} for event/no event.
#'
#' @param fun_aggregate Optional function used to aggregate multiple values
#'   within the same interval (e.g., \code{mean}, \code{median}). For
#'   \code{method = "measurement"}, this determines how measurements are
#'   summarized.
#'
#' @param fill Value used to fill missing cells in the wide-format output.
#'
#' @param ... Additional arguments passed to internal computations.
#'
#' @return A \code{data.table} in wide format with one row per subject and
#'   one column per interval. Column names are constructed as
#'   \code{paste0(name, "_", interval)}.
#'
#' @details
#' The function works by aligning long-format observations with a discrete
#' time grid using either rolling joins (for events and measurements) or
#' interval overlap calculations (for exposures).
#'
#' It is designed to be composable and can be used as a building block for
#' custom mapping functions supplied to \code{long_to_wide()}.
#'
#' @examples
#' library(data.table)
#'
#' ## Construct a simple discrete time grid
#' grid <- data.table(
#'   id = rep(1:2, each = 3),
#'   interval = rep(0:2, times = 2),
#'   start_interval = rep(c(0, 1, 2), times = 2),
#'   end_interval   = rep(c(1, 2, 3), times = 2)
#' )
#'
#' ## Example 1: event indicator carried forward
#' event_data <- data.table(
#'   id = c(1, 2),
#'   date = c(1.5, 0.5)
#' )
#'
#' discretize(
#'   method = "event",
#'   data = event_data,
#'   grid = grid,
#'   name = "Y",
#'   id = "id"
#' )
#'
#' ## Example 2: exposure percent within intervals
#' exposure_data <- data.table(
#'   id = c(1, 1, 2),
#'   start_date = c(0.2, 1.2, 0.0),
#'   end_date   = c(0.8, 2.5, 1.5)
#' )
#'
#' discretize(
#'   method = "exposure_percent",
#'   data = exposure_data,
#'   grid = grid,
#'   name = "A",
#'   id = "id"
#' )
#'
#'
#' @export
discretize <- function(method,
                       data,
                       grid,
                       name,
                       id,
                       threshold = NULL,
                       lookback_window = Inf,
                       values = c(1, 0),
                       fun_aggregate = NULL,
                       fill = NA,
                       ...) {
    value = exposure = start_date = end_date = start_followup_date = NULL
    x.interval = start_interval = end_interval = interval = date = NULL
    data <- copy(data)
    grid <- copy(grid)
    if (length(data) == 0) return(NULL)
    if (method == "measurement" && (length(fun_aggregate) == 0 || fun_aggregate != "last")){
        data[, interval := NA_integer_]
        # exact zero goes to interval 0
        data[date == 0, interval := 0L]
        # positive dates matched to (start_interval, end_interval]
        overlap <- data[date > 0,
                        interval := grid[data[date > 0],
                                         on = list(id,
                                                   start_interval < date,
                                                   end_interval >= date),
                                         x.interval]]
        ## discard information beyond the last interval (outside the grid)
        overlap <- overlap[!is.na(interval)]
        ## if (missing(fun_aggregate) || is.null(fun_aggregate)){
        if (is.function(fun_aggregate)){
            overlap[,list(value = fun_aggregate(value)),by = c(id,"interval")]
        }
    }
    if ((method %chin% c("event","locf","time_since_event"))
        || (method == "measurement" && (length(fun_aggregate)>0 &&fun_aggregate == "last"))){
        setnames(grid,"end_interval","date")
        if (method == "event" && length(values) == 2) {
            data[, value := values[[1]]]
        } else {
            if (method == "time_since_event"){
                set(data,j = "event_date",value = data[["date"]])
            }else{
                stopifnot("value" %in% names(data))
            }
        }
        # map observations/events to the grid using rolling join
        setkeyv(grid, c(id, "date"))
        setkeyv(data, c(id, "date"))
        overlap <- data[grid, roll = lookback_window]
        if (method == "time_since_event"){
            if (is.function(fun_aggregate)){
                set(overlap,j = "value",
                    value = fun_aggregate(as.numeric(overlap[["date"]]-overlap[["event_date"]])))
            }else{
                set(overlap,j = "value",
                    value = as.numeric(overlap[["date"]]-overlap[["event_date"]]))
            }
        }
        if (length(values) == 2) {
            # no event observed within roll window
            overlap[is.na(value), value := values[[2]]]
        }
    }
    # FIXME: in equidistant grids methods 'event' and 'any_exposure' can be done without
    #        foverlaps by using roll or even faster vector comparison of dates
    if (method %chin% "event_interval") {
        setnames(data, "date", "start_date")
        set(data, j = "end_date", value = data[["start_date"]])
    }
    if (method %chin%c("exposure_time","exposure_percent","event_interval","any_exposure","has_exposure")) {
        # calculate overlap of intervals in long data with grid intervals  
        setkeyv(data, c(id, "start_date", "end_date"))
        setkeyv(grid, c(id, "start_interval", "end_interval"))
        overlap <- data.table::foverlaps(x = grid,
                                         by.x = c(id, "start_interval", "end_interval"),
                                         y = data,
                                         by.y = c(id, "start_date", "end_date"),
                                         type = "any",
                                         nomatch = NA)
        if (method == "event_interval") {
            overlap[, exposure :=
                          1 * (
                              (interval == 0 & start_date == 0) |
                              (start_date <= end_interval & start_date >= start_interval)
                          )]
            overlap[is.na(exposure), exposure := 0]
            overlap <- overlap[, list(value = 1 * any(exposure == 1)), by = c(id, "interval")]
        } else {
            overlap[, exposure := (pmin(end_interval, end_date) -
                                   pmax(start_interval, start_date))]
            overlap[is.na(exposure), exposure := 0]
            # FIXME: baseline exposure is often only the start of exposure
            #        but in order to have adherence to the regimen we need to set it
            if ("start_followup_date" %chin%names(overlap)){
                overlap[interval == 0 & start_date == start_followup_date, exposure := 1]
            }
            setkeyv(overlap, c(id, "interval"))
            if (method %chin% c("exposure_time","exposure_percent","any_exposure","has_exposure")){
                if (method == "exposure_percent"){
                    # FIXME: the length of interval 0 is zero
                    overlap <- rbind(
                        overlap[interval == 0, list(interval = interval,value = sum(exposure)),by = id],
                        overlap[interval>0, list(value = sum(exposure)/(end_interval[1]-start_interval[1])), by = c(id, "interval")]
                    )
                }else{
                    overlap <- overlap[, list(value = sum(exposure)), by = c(id, "interval")]
                }
            }
            if (method %chin% c("has_exposure","any_exposure")){
                overlap[, value := 1 * (value > threshold)]
            }
        }
    }
    setkeyv(overlap, c(id, "interval"))
    wide <- fast_cast(x = overlap,
                      name = name,
                      id = id,
                      value_col = "value",
                      fill = fill,
                      fun_aggregate = fun_aggregate)        
    setkey(wide, id)
    setnames(wide, c(id, paste0(name, "_", names(wide)[-1])))
    setkeyv(wide, id)
    wide[]
}


######################################################################
### discretize_exposure.R ends here
