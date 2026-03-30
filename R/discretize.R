### discretize_exposure.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 24 2026 (11:04) 
## Version: 
## Last-Updated: mar 30 2026 (14:32) 
##           By: Thomas Alexander Gerds
##     Update #: 76
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
    value = exposure = start_date = end_date = NULL
    x.interval = start_interval = end_interval = interval = date = NULL
    data <- copy(data)
    grid <- copy(grid)
    if (length(data) == 0) return(NULL)
    if ((method == "measurement" && fun_aggregate != "last")){
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
        ## if (missing(fun_aggregate) || is.null(fun_aggregate)){
        if (is.function(fun_aggregate)){
            overlap[,list(value = fun_aggregate(value)),by = c(id,"interval")]
        }
    }
    if ((method %chin% c("event","locf"))
        || (method == "measurement" && fun_aggregate == "last")){
        setnames(grid,"end_interval","date")
        if (method == "event" && length(values) == 2) {
            data[, value := values[[1]]]
        } else {
            stopifnot("value" %in% names(data))
        }
        # map observations/events to the grid using rolling join
        setkeyv(grid, c(id, "date"))
        setkeyv(data, c(id, "date"))
        overlap <- data[grid, roll = lookback_window]
        if (length(values) == 2) {
            # no event observed within roll window
            overlap[is.na(value), value := values[[2]]]
        }
    }
    # FIXME: in equidistant grids methods 'event' and 'any_exposure' can be done without
    #        foverlaps by using roll or even faster vector comparison of dates
    if (method %chin% c("event","event_interval")) {
        setnames(data, "date", "start_date")
        set(data, j = "end_date", value = data[["start_date"]])
    }    
    if (method %in%c("exposure_time","event_interval","any_exposure")) {
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
            # baseline exposure
            overlap[interval == 0 & start_date == 0, exposure := 1]
            setkeyv(overlap, c(id, "interval"))
            overlap <- overlap[, list(value = sum(exposure)), by = c(id, "interval")]
            if (is.numeric(threshold)) {
                overlap[, value := 1 * (value > threshold)]
            }
        }
    }
    setkeyv(overlap, c(id, "interval"))
    wide <- fast_cast(overlap,
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
