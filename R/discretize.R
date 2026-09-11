### discretize.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 22 2024 (14:07)
## Version:
## Last-Updated: sep  8 2026 (18:07)
##           By: Thomas Alexander Gerds
##     Update #: 396
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
#' Discretize longitudinal data onto a discrete time grid
#'
#' Prepares long-format longitudinal data for analysis on a discrete time grid.
#' This is the high-level conversion function; \code{\link{map_data_to_grid}}
#' performs the individual data-to-grid mappings.
#'
#' \code{start_followup_date} is the calendar date at which follow-up starts, or
#' zero when the supplied dates are already on the time-on-study scale.
#'
#' @param x An object of class \code{"rtmle"} containing long-format data added
#'   with \code{\link{add_long_data}} and baseline data added with
#'   \code{\link{add_baseline_data}}.
#' @param start_followup_date Character scalar naming the Date-valued baseline
#'   column that contains each subject's follow-up start. That date is
#'   subtracted from outcome, competing-risk, censoring, and time-varying
#'   covariate dates. Use \code{0} (or omit the argument) when all long-format
#'   dates are already durations on the time-on-study scale. If long-format
#'   dates are \code{Date} objects, this argument must name a \code{Date} column;
#'   the long-format date columns are converted to numeric days since each
#'   subject's start.
#' @param mappings Named list of instructions for how to map long-format variables
#'   onto the discrete time scale. Names of \code{mappings} are output variable
#'   names. Each mapping may optionally contain \code{variable}, naming the
#'   source variable in \code{x$long_data$timevar_data}. If omitted, the source
#'   variable is assumed to have
#'   the same name as the output variable.
#'
#'   Built-in methods include \code{"measurement"}, \code{"locf"},
#'   \code{"event"}, \code{"event_interval"}, \code{"any_exposure"},
#'   \code{"has_exposure"}, \code{"exposure_time"}, and
#'   \code{"exposure_percent"}.
#'
#' @param verbose Logical. If \code{FALSE}, suppress informational messages
#'   generated while checking date formats.
#' @param ... alternative way to specify elements of \code{mappings}.
#' @details The function discretizes dates of events and concomitant marker
#' information. Multiple wide-format variables may result from a single
#' long-format variable. Calendar dates are converted to numeric time since
#' subject-specific follow-up start before the discrete grid is constructed.
#' The conversion is recorded in the object so that repeated calls do not
#' subtract the start date a second time; supplying new data with
#' \code{\link{add_long_data}} resets this state.
#' @return The modified \code{rtmle} object with wide-format outcome and
#'   time-varying covariate data stored in \code{x$data}.
#' @seealso \code{\link{rtmle_init}}, \code{\link{add_baseline_data}},
#'   \code{\link{add_long_data}}, \code{\link{map_data_to_grid}},
#'   \code{\link{long_to_wide}},
#'   \code{\link{prepare_rtmle_data}}
#' @examples
#' set.seed(17)
#' x <- rtmle_init(time_grid = 0:2,name_id = "id",
#'                 name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- add_long_data(x,outcome_data=data.frame(id=c(1,3,5),date=c(0.1,1.3,1.9)),
#'                      censored_data=data.frame(id=c(2,6),date=c(0.5,1.5)),
#'                     competing_data=data.frame(id=c(4,7),date=c(1.1,1.8)),
#'                     timevar_data=list(A=data.frame(id=c(1,1,1,2,2,3,4,4),
#'                                            start_date=c(0,.3,1.3,0,1,0,1.1,1.5),
#'                                              end_date=c(.3,.6,2,1,2,.5,1.4,1.9)),
#'                                       V=data.frame(id=c(1,1,1,2,2,3,4,4),
#'                                            start_date=c(0,.5,1,0,1.3,0,0,0.4),
#'                                              end_date=c(.3,.8,1.8,.5,1.4,2,.4,1.6)),
#'                                       L=data.frame(id=c(1,2,2,2,3,4,5,6,7,7),
#'                                                  date=c(0,0,.25,.75,0,0,0,0,0,1.4),
#'                                                 value=c(4,35,27.7,28.2,8.8,2,3.1,7,7.7,8.4))))
#' x <- add_baseline_data(x,data=data.frame(id=1:7,age=40:46))
#' x <- discretize(x,L=list(method="locf"),A=list(method="exposure_percent"))
#' x$data$timevar_data$L
#' x$data$timevar_data$V
#' x <- discretize(x,L=list(method="measurement",fun_aggregate="median"))
#' x$data
#' # multiple wide format variables from a single long format variable
#' x <- discretize(
#'                   x,
#'                   L=list(method="locf"),
#'                   A=list(method="exposure_percent"),
#'                   B=list(variable = "A", method="any_exposure"),
#'                   C=list(variable = "A", method="has_exposure",threshold=0.7),
#'                   start_followup_date=0
#' )
#' @export
discretize <- function(x,
                         start_followup_date,
                         mappings,
                         verbose = TRUE,
                         ...){
    start_interval = end_interval = interval = end_followup = censored_date =  competing_date = outcome_date = NULL
    id_column <- x$names$id
    breaks <- x$time_grid_scale

    if (length(x$data$baseline_data) == 0){
        stop("To discretize the long format data we need the baseline data stored as 'x$data$baseline_data' with the subject id variable.\nUse the function 'add_baseline_data' to add this.")
    }
    if (length(x$long_data) == 0) {
        stop("x$long_data has length 0.")
    }

    #
    # prepare mappings for all time-varying variables
    #
    Vnames <- names(x$long_data$timevar_data)
    if (length(Vnames) > 0){
        if (any(duplicated(Vnames))) {
            stop("Duplicated names found in names(x$long_data$timevar_data). Variables must have distinct names.")
        }
        if (missing(mappings)) mappings <- NULL
        dots <- list(...)
        if ((hit <- match("breaks", names(dots), nomatch = 0)) > 0){
            warning("rtmle::discretize: Argument 'breaks' is obsolete. Break points are stored in the object as x$time_grid and x$time_grid_scale.")
            dots <- dots[-hit]
        }
        if (length(dots) == 0) {
            dots <- NULL
        } else {
            bad <- is.null(names(dots)) || any(names(dots) == "")
            if (bad) stop("All ... arguments must be named.")
        }
        mappings <- c(mappings, dots)
        if (is.null(mappings)) {
            mappings <- list()
        }
        ## registry of built-in methods
        discretize_methods <- list(
            measurement      = list(method = "measurement",      fun = map_data_to_grid, columns = c("date","value"), lookback_window = Inf,fun_aggregate = "last"),
            locf             = list(method = "locf",             fun = map_data_to_grid, columns = c("date","value"), lookback_window = Inf),
            event            = list(method = "event",            fun = map_data_to_grid, columns = "date"),
            periodic_event   = list(method = "event_interval",   fun = map_data_to_grid, columns = "date"),
            time_since_event = list(method = "time_since_event", fun = map_data_to_grid, columns = "date"),
            chronic_disease  = list(method = "time_since_event", fun = map_data_to_grid, columns = "date", fun_aggregate = function(x){cut(x,breaks = c(-Inf,0,6*30.45,Inf),labels = c("never","acute","chronic"))}),
            event_interval   = list(method = "event_interval",   fun = map_data_to_grid, columns = "date"),
            any_exposure     = list(method = "any_exposure",     fun = map_data_to_grid, columns = c("start_date","end_date"), threshold = 0),
            has_exposure     = list(method = "has_exposure",     fun = map_data_to_grid, columns = c("start_date","end_date"), threshold = 0.5),
            exposure_time    = list(method = "exposure_time",    fun = map_data_to_grid, columns = c("start_date","end_date")),
            exposure_percent = list(method = "exposure_percent", fun = map_data_to_grid, columns = c("start_date","end_date"))
        )
        known_methods <- names(discretize_methods)

        # prepare mappings
        mappings <- lapply(names(mappings), function(Variable_name){
            spec <- mappings[[Variable_name]]
            if (is.character(spec) && length(spec) == 1L){
                spec <- list(method = spec)
            }
            if (!is.list(spec)) {
                stop("Mapping for '", Variable_name, "' must be a character string or a named list.")
            }
            long_format_varname <- spec$variable
            if (is.null(long_format_varname)) long_format_varname <- Variable_name
            if (!(long_format_varname %in% Vnames)){
                stop("Mapping '", Variable_name, "' refers to longformat variable '", long_format_varname,
                     "', which is not found in names(x$long_data$timevar_data).")
            }
            method_obj <- spec$method
            fun_obj <- spec$fun
            ## user can supply method as a function directly
            if (is.function(method_obj)) {
                fun_obj <- method_obj
                method_name <- Variable_name
            } else {
                method_name <- method_obj
            }
            if (is.character(method_name) && length(method_name) == 1L && method_name %in% known_methods){
                base <- discretize_methods[[method_name]]
                out <- c(spec, base)
                out <- out[!duplicated(names(out))]
                out$method <- method_name
                if (is.null(out$fun)) out$fun <- base$fun
                if (is.null(out$args)) out$args <- list()
            } else {
                ## user-defined mapping
                out <- spec
                if (is.null(fun_obj) || !is.function(fun_obj)){
                    stop("Unknown method for '", Variable_name,
                         "'. For custom mappings, supply either method=<function> or fun=<function>.")
                }
                out$fun <- fun_obj
                if (is.null(out$method)) out$method <- Variable_name
                if (is.null(out$columns)) {
                    stop("Custom mapping for '", Variable_name,
                         "' must specify required input columns via columns = c(...).")
                }
                if (is.null(out$args)) out$args <- list()
            }
            out$target <- Variable_name
            out$variable <- long_format_varname
            out
        })
        names(mappings) <- vapply(mappings, `[[`, "", "target")
        #
        ## validate required columns for each mapping
        #
        for (Variable_name in names(mappings)){
            m <- mappings[[Variable_name]]
            long_format_varname <- m$variable
            if (!(all(m$columns %chin% names(x$long_data$timevar_data[[long_format_varname]])))){
                stop(
                    paste0(
                        "To apply the long-to-wide mapping method '", m$method,
                        "' for output variable '", Variable_name,
                        "' using source variable '", long_format_varname, "',\n",
                        "the source data must contain the following columns:\n",
                        paste0(m$columns, collapse = ", "), "."
                    )
                )
            }
        }
        #
        ## auto-map source variables that were not used at all
        #
        mapped_sources <- unique(vapply(mappings, `[[`, "", "variable"))
        if (length(unmapped <- setdiff(Vnames, mapped_sources)) > 0){
            guessed <- lapply(unmapped, function(v){
                guessed_mapping <- NULL
                for (m in discretize_methods) {
                    if (all(m$columns %in% names(x$long_data$timevar_data[[v]]))) {
                        guessed_mapping <- m
                        break
                    }
                }
                if (is.null(guessed_mapping)) {
                    stop("Could not determine method for time-varying variable ", v,
                         " from the columns in 'data'; please specify the method explicitly.")
                }
                guessed_mapping$variable <- v
                guessed_mapping$target <- v
                guessed_mapping
            })
            names(guessed) <- unmapped
            mappings <- c(mappings, guessed)
        }
        # timevar date formats
        tv_date_formats <- vapply(Vnames, function(v){
            current_data <- x$long_data$timevar_data[[v]]
            if ("date" %in% names(current_data)){
                if (inherits(current_data[["date"]], "Date")) {
                    "Date"
                } else if (is.numeric(current_data[["date"]])) {
                    "numeric"
                } else {
                    class(current_data[["date"]])[[1L]]
                }
            } else {
                if (!all(c("start_date", "end_date") %in% names(current_data))) {
                    stop("Time-varying variable '", v,
                         "' must contain date or start_date/end_date columns.")
                }
                start_class <- if (inherits(current_data[["start_date"]], "Date")) {
                    "Date"
                } else if (is.numeric(current_data[["start_date"]])) {
                    "numeric"
                } else {
                    class(current_data[["start_date"]])[[1L]]
                }
                end_class <- if (inherits(current_data[["end_date"]], "Date")) {
                    "Date"
                } else if (is.numeric(current_data[["end_date"]])) {
                    "numeric"
                } else {
                    class(current_data[["end_date"]])[[1L]]
                }
                if (!identical(start_class, end_class)) {
                    stop("The start_date and end_date columns for time-varying variable '",
                         v, "' must have the same class.")
                }
                start_class
            }
        }, character(1L))

    } else {
        tv_date_formats <- NULL
    }
    # outcome date formats
    ocd_date_formats <- vapply(
        intersect(c("outcome_data", "censored_data", "competing_data"),
                  names(x$long_data)),
        function(v) {
            current_data <- x$long_data[[v]]
            if (!("date" %in% names(current_data))) {
                stop("The long-format ", v, " data must contain a 'date' column.")
            }
            if (inherits(current_data[["date"]], "Date")) {
                "Date"
            } else if (is.numeric(current_data[["date"]])) {
                "numeric"
            } else {
                class(current_data[["date"]])[[1L]]
            }
        },
        character(1L)
    )
    #
    # check if all date formats are the same
    #
    date_formats <- c(tv_date_formats, ocd_date_formats)
    if (length(unique(date_formats)) > 1){
        ctab <- data.table(
            "variable" = names(date_formats),
            "storage" = rep(c("long_data$timevar_data", "long_data"),
                            c(length(tv_date_formats), length(ocd_date_formats))),
            "class_of_date_variable" = unname(date_formats))
        if (isTRUE(verbose)) {
            cat("\n\nAnalysis of date variables:\n\n")
            print(ctab)
            cat("\n")
        }
        stop("All date variables in long format data need to have the same class (either numeric or Date).")
    }

    date_kind <- if (length(date_formats) == 0L) NULL else unname(date_formats[[1L]])
    if (!is.null(date_kind) && !date_kind %in% c("numeric", "Date")) {
        stop("Date variables in long format data must be numeric or Date.")
    }

    if (is.null(x$progress)) x$progress <- list()
    if (is.null(x$progress$substracted_start_followup_date)) {
        x$progress$substracted_start_followup_date <- FALSE
    }

    # A calendar-date history needs a start date for each subject. Convert all
    # event, measurement, and exposure dates to numeric time since that start
    # before constructing the grid. The flag makes repeated calls idempotent;
    # add_long_data() resets it when new long-format data are supplied.
    has_start_followup <- !missing(start_followup_date) &&
        length(start_followup_date) > 0L &&
        !(is.numeric(start_followup_date) &&
          length(start_followup_date) == 1L &&
          !is.na(start_followup_date[[1L]]) &&
          start_followup_date[[1L]] == 0)

    if (!has_start_followup) {
        if (identical(date_kind, "Date") &&
            !isTRUE(x$progress$substracted_start_followup_date)) {
            stop("Calendar-date long-format data require `start_followup_date` to name a Date column in the baseline data.")
        }
        if (missing(start_followup_date)) {
            x$diagnostics$missing_start_followup_variable <- "Assume 0 is time zero."
        }
        pop <- x$data$baseline_data[
            , c(id_column),
            with = FALSE
        ][, start_followup_date := rep(0, .N)]
    } else {
        if (length(start_followup_date) != 1L ||
            !is.character(start_followup_date) ||
            is.na(start_followup_date) ||
            !nzchar(start_followup_date) ||
            match(start_followup_date, names(x$data$baseline_data), nomatch = 0) == 0) {
            stop("Argument start_followup_date must be the name (as character) of a variable in x$data$baseline_data")
        }
        x$names$start_followup_date <- start_followup_date
        pop <- x$data$baseline_data[
            , c(id_column, start_followup_date),
            with = FALSE
        ]
        setnames(pop, start_followup_date, "start_followup_date")
        x$diagnostics$missing_start_followup <- NULL
    }

    if (identical(date_kind, "Date") &&
        !isTRUE(x$progress$substracted_start_followup_date)) {
        if (!has_start_followup ||
            !inherits(pop[["start_followup_date"]], "Date")) {
            stop("Calendar-date long-format data require a Date-valued `start_followup_date` column.")
        }
        start_dates <- pop[["start_followup_date"]]
        if (anyNA(start_dates)) {
            stop("The start_followup_date column cannot contain missing values.")
        }
        subtract_start_followup <- function(current_data, date_columns, label) {
            if (is.null(current_data) || NROW(current_data) == 0L) {
                return(invisible(NULL))
            }
            row_index <- match(current_data[[id_column]], pop[[id_column]])
            if (anyNA(row_index)) {
                stop("The ", label,
                     " data contain subject ids that are not present in the baseline data.")
            }
            for (date_column in date_columns) {
                if (date_column %in% names(current_data)) {
                    data.table::set(
                        current_data,
                        j = date_column,
                        value = as.numeric(
                            current_data[[date_column]] - start_dates[row_index]
                        )
                    )
                }
            }
            invisible(NULL)
        }
        for (data_name in intersect(
            c("outcome_data", "censored_data", "competing_data"),
            names(x$long_data)
        )) {
            subtract_start_followup(
                x$long_data[[data_name]],
                "date",
                paste0("long_data$", data_name)
            )
        }
        for (data_name in names(x$long_data$timevar_data)) {
            subtract_start_followup(
                x$long_data$timevar_data[[data_name]],
                c("date", "start_date", "end_date"),
                paste0("long_data$timevar_data$", data_name)
            )
        }
        pop[, start_followup_date := rep(0, .N)]
        x$progress$substracted_start_followup_date <- TRUE
    } else if (isTRUE(x$progress$substracted_start_followup_date)) {
        if (!is.null(date_kind) && !identical(date_kind, "numeric")) {
            stop("Long-format dates marked as converted must be numeric.")
        }
        pop[, start_followup_date := rep(0, .N)]
    } else if (has_start_followup &&
               inherits(pop[["start_followup_date"]], "Date")) {
        stop("A Date-valued start_followup_date requires Date-valued long-format dates.")
    }

    # Keep the source long-format tables in their original column naming while
    # constructing the follow-up history from private copies.
    outcome_data <- NULL
    if ("outcome_data" %in% names(x$long_data)) {
        outcome_data <- data.table::copy(x$long_data$outcome_data)
        if (!"date" %in% names(outcome_data)) {
            stop("The long-format outcome_data must contain a 'date' column.")
        }
        setnames(outcome_data, "date", "outcome_date")
    } else {
        stop("discretize: Object does not contain outcome data.")
    }
    competing_data <- NULL
    if (length(x$names$competing) > 0 &&
        "competing_data" %in% names(x$long_data) &&
        NROW(x$long_data$competing_data) > 0L) {
        competing_data <- data.table::copy(x$long_data$competing_data)
        setnames(competing_data, "date", "competing_date")
    }
    censored_data <- NULL
    if (length(x$names$censoring) > 0 &&
        "censored_data" %in% names(x$long_data) &&
        NROW(x$long_data$censored_data) > 0L) {
        censored_data <- data.table::copy(x$long_data$censored_data)
        setnames(censored_data, "date", "censored_date")
    }

    # Calculate each subject's observed end of follow-up.
    if (!is.null(competing_data)) {
        pop <- competing_data[pop, on = id_column]
    } else {
        pop[, competing_date := Inf]
    }
    if (!is.null(censored_data)) {
        pop <- censored_data[pop, on = id_column]
    } else {
        pop[, censored_date := Inf]
    }
    pop <- outcome_data[pop, on = id_column]
    pop[, end_followup := pmin(
        censored_date,
        competing_date,
        outcome_date,
        na.rm = TRUE
    )]
    if (any(is.na(pop$end_followup))) {
        stop("Missing values in end of followup information")
    }
    if (NROW(pop) > 0L &&
        max(pop$end_followup) < rev(x$time_grid_scale)[2]) {
        stop(paste0(
            "The maximal followup time in the data is ",
            max(pop$end_followup), ",\n",
            "but the last interval on the time grid starts at ",
            rev(x$time_grid_scale)[2], ".\n",
            "Please remove all time-grid values that are beyond the maximal followup time."
        ))
    }

    # An event after a competing event or censoring is not observed, and a
    # competing/censoring event after the outcome is irrelevant. Remove these
    # dates before mapping the outcome history to the wide representation.
    pop[
        !is.na(outcome_date) & outcome_date > end_followup,
        outcome_date := NA_real_
    ]
    pop[
        !is.na(censored_date) & !is.na(competing_date) &
            competing_date > censored_date,
        competing_date := NA_real_
    ]
    pop[
        !is.na(competing_date) & !is.na(outcome_date) &
            competing_date >= outcome_date,
        competing_date := NA_real_
    ]
    pop[
        !is.na(censored_date) & !is.na(outcome_date) &
            censored_date >= outcome_date,
        censored_date := NA_real_
    ]
    pop[
        !is.na(censored_date) & !is.na(competing_date) &
            censored_date >= competing_date,
        censored_date := NA_real_
    ]

    # Build the subject-specific discrete grid on the numeric time-on-study
    # scale.
    grid <- pop[
        , list(
            start_followup_date = start_followup_date,
            end_interval = start_followup_date + breaks,
            end_followup = end_followup
        ),
        by = id_column
    ]
    grid[
        , start_interval := c(start_followup_date[[1L]], end_interval[-.N]),
        by = id_column
    ]
    grid <- grid[start_interval <= end_followup]
    grid[, end_followup := NULL]
    grid[, interval := seq_len(.N) - 1L, by = id_column]
    data.table::setcolorder(
        grid,
        c(id_column, "interval", "start_interval", "end_interval")
    )

    outcome_for_widen <- pop[
        !is.na(outcome_date) & is.finite(outcome_date),
        c(id_column, "outcome_date"),
        with = FALSE
    ]
    if (NROW(outcome_for_widen) == 0L) {
        stop("discretize: Object does not contain outcome dates or all outcome dates are administratively censored.")
    }
    competing_for_widen <- if (!is.null(competing_data)) {
        pop[!is.na(competing_date) & is.finite(competing_date),
            c(id_column, "competing_date"), with = FALSE]
    } else {
        NULL
    }
    censored_for_widen <- if (!is.null(censored_data)) {
        pop[!is.na(censored_date) & is.finite(censored_date),
            c(id_column, "censored_date"), with = FALSE]
    } else {
        NULL
    }

    # map outcome data
    x$data$outcome_data <- widen_outcome(
        x,
        outcome_data = outcome_for_widen,
        censored_data = censored_for_widen,
        competing_data = competing_for_widen,
        grid = grid,
        fun_aggregate = NULL
    )
    #
    # map timevarying data
    #
    if (length(Vnames) > 0){
        for (Variable_name in names(mappings)){
            m <- mappings[[Variable_name]]
            long_format_varname <- m$variable
            fun <- m$fun
            m[["fun"]] <- m[["variable"]] <- NULL
            if(m$method %chin% c("event","locf") && length(m$lookback_window) == 0){
                lookback_window <- Inf
            }
            args <- c(m,
                      list(
                          data = x$long_data$timevar_data[[long_format_varname]],
                          grid = grid,
                          name = Variable_name,
                          id = id_column,
                          values = c(1, 0),
                          fill = NA
                      ))
            args <- args[!duplicated(names(args))]
            ## keep only args accepted by fun unless it has ...
            fml <- names(formals(fun))
            if (!("..." %chin% fml)){
                args <- args[intersect(names(args), fml)]
            }
            x$data$timevar_data[[Variable_name]] <- do.call(fun, args)
            ## TODO: Uncomment when OK
            ## ## Add mappings names for each constructed variable to x
            ## x$discretize_mappings[[Variable_name]] <- mappings[[Variable_name]]$method
        }
    }
    x
}

######################################################################
### discretize.R ends here
