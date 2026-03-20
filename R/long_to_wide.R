### long_to_wide.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 22 2024 (14:07) 
## Version: 
## Last-Updated: mar 14 2026 (07:02) 
##           By: Thomas Alexander Gerds
##     Update #: 197
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# Argument 'fun': if more than one instance (such as treatment) happens
#        within one interval
#        an aggregate function determines the value inside
#        the interval
#' Transform data from long to wide format
#'
#' This function is used to prepare a discretized time analysis of longituginal data.
#' 
#' The start_followup_date is the calendar date where the follow up
#' starts or can be zero if the followup data are readily on the time
#' on study scale
#' @param x object of class \code{rtmle}
#' @param start_followup_date Start of followup date which is substracted from the dates of
#' the outcomes, competing risks, censoring (end-of-followup), and time-varying covariates.
#' If missing it is assumed to be zero for all in which case the dates of the outcomes,
#' competing risks, censoring (end-of-followup), and time-varying covariates must be
#' given on the time on study scale.
#' @param fun function used to map time-dependent covariate information onto the intervals.
#' Can be a single function or a named list with 
#' @param verbose Logical. If \code{FALSE} suppress all messages. 
#' @details The function discretizes data on a given time grid. Dates of events in long format. 
#' @examples
#' set.seed(17)
#' tau <- 3
#' ld <- simulate_long_data(n = 91,number_visits = 20,
#'                          beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),
#'                          register_format = TRUE)
#' 
#' x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- add_long_data(x,outcome_data=ld$outcome_data,censored_data=ld$censored_data,
#'                     competing_data=ld$competing_data,
#'                     timevar_data=ld$timevar_data)
#' x <- add_baseline_data(x,data=ld$baseline_data)
#' x <- long_to_wide(x,breaks = seq(0,2000,30.45*12))
#'                   
#' x
#' @export
long_to_wide <- function(x,
                         start_followup_date,
                         fun = function(x){1*(sum(x)>0)},
                         verbose = TRUE){
    previous_date = interval = end_followup = censored_date =  competing_date = outcome_date = NULL
    breaks = x$time_grid
    if (length(x$long_data) == 0) {return(NULL)}
    Vnames <- names(x$long_data$timevar_data)
    if (any(duplicated(Vnames))) stop("Duplicated names found in names(x$long_data$timevar_data). Variables must have distinct names.")
    # check if all dates have the same format
    tv_date_formats <- sapply(Vnames,function(v){
        if ("date"%in%names(x$long_data$timevar_data[[v]])){
            class(x$long_data$timevar_data[[v]][["date"]])
        }else{
            uclass <- unique(c(class(x$long_data$timevar_data[[v]][["start_exposure"]]),
                               class(x$long_data$timevar_data[[v]][["end_exposure"]])))
            stopifnot(length(uclass) == 1)
            uclass
        }
    })
    ocd_date_formats <- sapply(intersect(c("outcome_data","censored_data","competing_data"),names(x$long_data)),function(v){
        class(x$long_data[[v]][["date"]])
    })
    if (length(unique(c(tv_date_formats,ocd_date_formats))) > 1){
        ctab <- data.table(
            "variable" = names(c(tv_date_formats,ocd_date_formats)),
            "storage" = rep(c("long_data$timevar_data","long_data"),c(length(tv_date_formats),length(ocd_date_formats))),
            "class_of_date_variable" = c(tv_date_formats,ocd_date_formats))
        cat("\n\nAnalysis of date variables:\n\n")
        print(ctab)
        cat("\n")
        stop(paste0("All date variables in long format data need to have the same class (either numeric or Date)."))
    }
    if (length(x$data$baseline_data) == 0){
        stop("To discretize the long format data we need the baseline data stored as 'x$data$baseline_data' with the subject id variable.\nUse the function 'add_baseline_data' to add this.")
    }
    if (missing(start_followup_date)||start_followup_date[1] == 0){
        if (missing(start_followup_date)){
            x$details <- c(x$details,list("Missing start_followup_date variable, for now assume 0, which implies that all event times must be given on the scale: 'time since 0'."))
        }
        pop <- x$data$baseline_data[,c(x$names$id),with = FALSE][,start_followup_date := rep(0,.N)]
    }else{
        if (!is.character(start_followup_date) || match(start_followup_date,names(x$data$baseline_data),nomatch = 0) == 0){
            stop("Argument start_followup_date must be the name (as character) of a variable in x$data$baseline_data")
        }
        x$names$start_followup_date <- start_followup_date
        pop <- x$data$baseline_data[,c(x$names$id,start_followup_date),with = FALSE]
        setnames(pop,start_followup_date,"start_followup_date")
    }
    #
    # outcome, censored and competing risk define end-of-followup
    #
    if (length(x$names$competing)>0 & length(x$long_data$competing_data)>0){
        if (length(x$long_data$outcome_data)>0){
            if (!("competing_date" %in% names(x$long_data$competing_data)))
                setnames(x$long_data$competing_data,"date","competing_date")
        }
        pop <- x$long_data$competing_data[pop,on = x$names$id]
    }else{
        pop[,competing_date := as.Date(Inf)]
    }
    if (length(x$names$censoring)>0 & length(x$long_data$censored_data)>0){
        if (!("censored_date" %in% names(x$long_data$censored_data)))
            setnames(x$long_data$censored_data,"date","censored_date")
        pop <- x$long_data$censored_data[pop,on = x$names$id]
    }else{
        pop[,censored_date := as.Date(Inf)]
    }
    if (!("outcome_date" %in% names(x$long_data$outcome_data)))
        setnames(x$long_data$outcome_data,"date","outcome_date")
    pop <- x$long_data$outcome_data[pop,on = x$names$id]
    pop[,end_followup := pmin(censored_date,competing_date,outcome_date,na.rm = TRUE)]
    pop[,c(x$names$id,"start_followup_date","end_followup"),with = FALSE]
    if (any(is.na(pop$end_followup)))stop("Missing values in end of followup information")
    id_colname <- x$names$id
    grid <- pop[,data.table::data.table(date=start_followup_date+breaks, end_followup = end_followup),by=id_colname]
    # FIXME: check for data after the end of followup?
    grid[,previous_date := c(0,date[-.N]),by = id_colname]
    grid <- grid[previous_date<=end_followup]
    grid[,end_followup := NULL]
    grid[,interval:=0:(.N-1),by=id_colname]
    grid <- pop[,.SD,.SDcols = c(x$names$id)][grid,on = id_colname]
    # now awkwardly reset the names
    if (length(x$names$competing)>0 & length(x$long_data$competing_data)>0){
        setnames(x$long_data$competing_data,"competing_date","date")
    }
    if (length(x$names$censoring)>0 & length(x$long_data$censored_data)>0){
        setnames(x$long_data$censored_data,"censored_date","date")
    }
    setnames(x$long_data$outcome_data,"outcome_date","date")
    # mapping outcome information to discrete time scale
    x$data$outcome_data <- widen_outcome(x,grid = grid,fun_aggregate = NULL)
    # time dependent variables including treatment
    for (Vname in Vnames){
        if (is.list(fun)){
            if (Vname %in% names(fun)){
                # use variable specific function
                vfun <- fun[[Vname]]
            } else{
                # default to last value carried forward for non-binary variables
                if ("value" %in% names(x$long_data$timevar_data[[Vname]])){
                    vfun <- function(x){x}
                }else{
                    # default to has positive value for non-binary variables
                    # note: in this case the function map_grid needs the values sorted 
                    vfun <- function(x){1*(sum(x)>0)}
                }
            }
        } else {
            vfun <- fun
        }
        stopifnot(is.function(vfun))
        # Deal with 3 different formats of the data
        # case 1: id, date (exposed = there is a date in the interval)
        # case 2: id, date, value (last value carried forward)
        # case 3: id, start, end (calculate overlap with interval)
        if ("value" %in% names(x$long_data$timevar_data[[Vname]])){
            # case 2: last value carried forward
            x$data$timevar_data[[Vname]] <- map_grid(grid=grid,
                                                     data=x$long_data$timevar_data[[Vname]],
                                                     name=Vname,
                                                     fun_aggregate = vfun,
                                                     values = NULL, 
                                                     rollforward=Inf,
                                                     id = x$names$id)
        } else{
            if (all(c("start_exposure","end_exposure") %in% names(x$long_data$timevar_data[[Vname]]))){
                # case 3: calculate overlap
                gg = copy(grid)
                setnames(gg,c("date","previous_date"),c("end_interval","start_interval"))
                x$data$timevar_data[[Vname]] <- discretize_timevarying_exposure(data = x$long_data$timevar_data[[Vname]],
                                                                                grid = gg,
                                                                                name = Vname,
                                                                                point_exposure = FALSE,
                                                                                threshold = 0,
                                                                                id = x$names$id)
            }else{ # case 1
                gg = copy(grid)
                setnames(gg,c("date","previous_date"),c("end_interval","start_interval"))
                x$data$timevar_data[[Vname]] <- discretize_timevarying_exposure(data = x$long_data$timevar_data[[Vname]],grid = gg,name = Vname,point_exposure = TRUE,id = x$names$id)
                # FIXME: when intervals are not equidistant the following does not work.
                ## length_interval=unique(round(diff(breaks),0))
                ## x$data$timevar_data[[Vname]] <- map_grid(grid=grid,data=x$long_data$timevar_data[[Vname]],name=Vname,values = c(1,0),fun_aggregate = vfun,rollforward=length_interval,id = x$names$id)
            }
        }
    }
    return(x)
}



######################################################################
### long_to_wide.R ends here
