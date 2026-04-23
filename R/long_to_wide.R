### long_to_wide.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 22 2024 (14:07) 
## Version: 
## Last-Updated: apr 23 2026 (17:40) 
##           By: Thomas Alexander Gerds
##     Update #: 388
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Transform data from long to wide format
#'
#' This function is used to prepare a discretized time analysis of longituginal data.
#' 
#' The start_followup_date is the calendar date where the follow up
#' starts or can be zero if the followup data are readily on the time
#' on study scale
#' @param x object of class \code{rtmle}
#' @param start_followup_date Date where the followup starts. This date is substracted from the dates of
#' the outcomes, competing risks, censoring (end-of-followup), and dates where the time-varying covariates change.
#' If missing it is assumed to be zero for all subjects in which case the dates of the outcomes,
#' competing risks, censoring (end-of-followup), and time-varying covariates must be
#' given on the time on study scale.
#' @param mappings Named list of instructions for how to map long-format variables
#'   onto the discrete time scale. Names of `mappings` are output variable names.
#'   Each mapping may optionally contain `variable`, naming the source variable in
#'   `x$long_data$timevar_data`. If omitted, the source variable is assumed to have
#'   the same name as the output variable.
#'
#'   Built-in methods include "measurement", "locf", "event", "event_interval",
#'   "any_exposure", "has_exposure", "exposure_time", and "exposure_percent".
#'
#' @param verbose Logical. If \code{FALSE} suppress all messages.
#' @param ... alternative way to specify elements of \code{mappings}.
#' @details The function discretizes dates of events and concomittant marker information.
#' Multiple wide format variables may result from a single long-format variable.
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
#' x <- long_to_wide(x,L=list(method="locf"),A=list(method="exposure_percent"))
#' x$data$timevar_data$L
#' x$data$timevar_data$V
#' x <- long_to_wide(x,L=list(method="measurement",fun_aggregate="median"))
#' x$data
#' # multiple wide format variables from a single long format variable
#' x <- long_to_wide(
#'                   x,
#'                   L=list(method="locf"),
#'                   A=list(method="exposure_percent"),
#'                   B=list(variable = "A", method="any_exposure"),
#'                   C=list(variable = "A", method="has_exposure",threshold=0.7),
#'                   start_followup_date=0
#' )
#' @export
long_to_wide <- function(x,
                         start_followup_date,
                         mappings,
                         verbose = TRUE,
                         ...){
    start_interval = end_interval = interval = end_followup = censored_date =  competing_date = outcome_date = NULL
    id_column <- x$names$id
    breaks <- x$time_grid_labels

    if (length(x$data$baseline_data) == 0){
        stop("To discretize the long format data we need the baseline data stored as 'x$data$baseline_data' with the subject id variable.\nUse the function 'add_baseline_data' to add this.")
    }
    if (length(x$long_data) == 0) {
        stop("x$long_data has length 0.")
    }

    #
    # start of followup
    #
    if (missing(start_followup_date) || start_followup_date[1] == 0){
        # FIXME: want to store this somewhere or simply print it as a message
        if (missing(start_followup_date)){
            message("Missing start_followup_date variable, for now assume 0.\nThis makes only sense if all event times are given on the scale: 'time since 0'.")
        }
        pop <- x$data$baseline_data[, c(id_column), with = FALSE][, start_followup_date := rep(0, .N)]
    } else {
        if (!is.character(start_followup_date) || match(start_followup_date, names(x$data$baseline_data), nomatch = 0) == 0){
            stop("Argument start_followup_date must be the name (as character) of a variable in x$data$baseline_data")
        }
        x$names$start_followup_date <- start_followup_date
        pop <- x$data$baseline_data[, c(id_column, start_followup_date), with = FALSE]
        setnames(pop, start_followup_date, "start_followup_date")
        x$diagnostics$missing_start_followup <- NULL
    }
    #
    # rename the date columns in order to calculate the end-of-followup
    #
    if (length(x$names$competing) > 0 && length(x$long_data$competing_data) > 0){
        setnames(x$long_data$competing_data, "date", "competing_date")
        on.exit(try(setnames(x$long_data$competing_data, "competing_date", "date"),silent = TRUE))
    }
    if (length(x$names$censoring) > 0 && length(x$long_data$censored_data) > 0){
        setnames(x$long_data$censored_data, "date", "censored_date")
        on.exit(try(setnames(x$long_data$censored_data, "censored_date", "date"),silent = TRUE))
    }
    setnames(x$long_data$outcome_data, "date", "outcome_date")
    on.exit(try(setnames(x$long_data$outcome_data, "outcome_date", "date"),silent = TRUE))
    #
    # competing risks
    # 
    if (length(x$names$competing) > 0 && length(x$long_data$competing_data) > 0){
        pop <- x$long_data$competing_data[pop, on = id_column]
    } else {
        pop[, competing_date := as.Date(Inf)]
    }
    #
    # right censoring
    # 
    if (length(x$names$censoring) > 0 && length(x$long_data$censored_data) > 0){
        pop <- x$long_data$censored_data[pop, on = id_column]
    } else {
        pop[, censored_date := as.Date(Inf)]
    }
    #
    # outcome
    #
    pop <- x$long_data$outcome_data[pop, on = id_column]
    pop[, end_followup := pmin(censored_date, competing_date, outcome_date, na.rm = TRUE)]
    pop[, c(id_column, "start_followup_date", "end_followup"), with = FALSE]
    if (any(is.na(pop$end_followup))) stop("Missing values in end of followup information")
    #
    # discrete time grid
    #
    grid <- pop[, list(start_followup_date = start_followup_date,
                       end_interval = start_followup_date + breaks,
                       end_followup = end_followup), by = id_column]
    grid[, start_interval := c(start_followup_date[1],end_interval[-.N]), by = id_column]
    grid <- grid[start_interval <= end_followup]
    grid[, end_followup := NULL]
    grid[, interval := 0:(.N - 1), by = id_column]
    grid <- pop[, .SD, .SDcols = c(id_column)][grid, on = id_column]
    data.table::setcolorder(grid, c(id_column, "interval", "start_interval", "end_interval"))
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
            warning("rtmle::long_to_wide: Argument 'breaks' is obsolete. Break points are stored in the object as x$time_grid and x$time_grid_labels.")
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
        long_to_wide_methods <- list(
            measurement      = list(method = "measurement",      fun = discretize, columns = c("date","value"), lookback_window = Inf,fun_aggregate = "last"),
            locf             = list(method = "locf",             fun = discretize, columns = c("date","value"), lookback_window = Inf),
            event            = list(method = "event",            fun = discretize, columns = "date"),
            periodic_event   = list(method = "event_interval",   fun = discretize, columns = "date"),
            time_since_event = list(method = "time_since_event", fun = discretize, columns = "date"),
            chronic_disease  = list(method = "time_since_event", fun = discretize, columns = "date", fun_aggregate = function(x){cut(x,breaks = c(-Inf,0,6*30.45,Inf),labels = c("never","acute","chronic"))}),
            event_interval   = list(method = "event_interval",   fun = discretize, columns = "date"),
            any_exposure     = list(method = "any_exposure",     fun = discretize, columns = c("start_date","end_date"), threshold = 0),
            has_exposure     = list(method = "has_exposure",     fun = discretize, columns = c("start_date","end_date"), threshold = 0.5),
            exposure_time    = list(method = "exposure_time",    fun = discretize, columns = c("start_date","end_date")),
            exposure_percent = list(method = "exposure_percent", fun = discretize, columns = c("start_date","end_date"))
        )
        known_methods <- names(long_to_wide_methods)
        
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
                base <- long_to_wide_methods[[method_name]]
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
                for (m in long_to_wide_methods) {
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
        tv_date_formats <- sapply(Vnames, function(v){
            if ("date" %in% names(x$long_data$timevar_data[[v]])){
                class(x$long_data$timevar_data[[v]][["date"]])
            } else {
                uclass <- unique(c(class(x$long_data$timevar_data[[v]][["start_date"]]),
                                   class(x$long_data$timevar_data[[v]][["end_date"]])))
                stopifnot(length(uclass) == 1)
                uclass
            }
        })

    } else {
        tv_date_formats <- NULL
    }
    # outcome date formats
    ocd_date_formats <- sapply(intersect(c("outcome_data","censored_data","competing_data"), names(x$long_data)), function(v){
        # here we rename the date variable back to date
        setnames(x$long_data[[v]],sub("data","date",v),"date")
        class(x$long_data[[v]][["date"]])
    })
    #
    # check if all date formats are the same
    #
    if (length(unique(c(tv_date_formats, ocd_date_formats))) > 1){
        ctab <- data.table(
            "variable" = names(c(tv_date_formats, ocd_date_formats)),
            "storage" = rep(c("long_data$timevar_data","long_data"), c(length(tv_date_formats), length(ocd_date_formats))),
            "class_of_date_variable" = c(tv_date_formats, ocd_date_formats))
        cat("\n\nAnalysis of date variables:\n\n")
        print(ctab)
        cat("\n")
        stop("All date variables in long format data need to have the same class (either numeric or Date).")
    }
    #
    # map outcome data
    #
    x$data$outcome_data <- widen_outcome(x, grid = grid, fun_aggregate = NULL)
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
        }
    }
    x
}

######################################################################
### long_to_wide.R ends here
