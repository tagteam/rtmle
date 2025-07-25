### long_to_wide.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 22 2024 (14:07) 
## Version: 
## Last-Updated: Jul 25 2025 (10:27) 
##           By: Thomas Alexander Gerds
##     Update #: 93
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
#' @param intervals a vector of time points that discretize the
#'     followup time into intervals
#' @param start_followup_date Start of followup date which is substracted from the dates of
#' the outcomes, competing risks, censoring (end-of-followup), and time-varying covariates.
#' If missing it is assumed to be zero for all in which case the dates of the outcomes,
#' competing risks, censoring (end-of-followup), and time-varying covariates must be
#' given on the time on study scale.
#' @param fun function used to map information onto the intervals. see Details.
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
#' x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
#'                   
#' x
#' @export
long_to_wide <- function(x,
                         intervals,
                         start_followup_date,
                         fun = function(x){1*(sum(x)>0)},
                         verbose = TRUE){
    interval = end_followup = censored_date =  competing_date = outcome_date = NULL
    if (length(x$long_data) == 0) {return(NULL)}
    Vnames <- names(x$long_data$timevar_data)
    if (any(duplicated(Vnames))) stop("Duplicated names found in names(x$long_data$timevar_data). Variables must have distinct names.")
    # check if all dates have the same format
    tv_date_formats <- sapply(Vnames,function(v){class(x$long_data$timevar_data[[v]][["date"]])})
    ocd_date_formats <- sapply(intersect(c("outcome_data","censored_data","competing_data"),names(x$long_data)),function(v){class(x$long_data[[v]][["date"]])})
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
    if (missing(intervals)) {stop("Need time intervals to map the long data.")}
    if (length(x$data$baseline_data) == 0){
        stop("To discretize the long format data we need the baseline data stored as 'x$data$baseline_data' with the subject id variable.\nUse the function 'add_baseline_data' to add this.")
    }
    if (missing(start_followup_date)){
        if (verbose) message("Missing start_followup_date variable, for now assume 0, which implies that all event times must be given on the time on study scale.")
        pop = x$data$baseline_data[,c(x$names$id),with = FALSE][,start_followup_date := rep(0,.N)]
    }else{
        if (!is.character(start_followup_date) || match(start_followup_date,names(x$data$baseline_data),nomatch = 0) == 0){
            stop("Argument start_followup_date must be the name (as character) of a variable in x$data$baseline_data")
        }
        pop <- x$data$baseline_data[,c(x$names$id,start_followup_date),with = FALSE]
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
    if (any(is.na(pop$end_followup)))stop("Missing values in end of followup information")
    grid <- pop[,data.table::data.table(date=start_followup_date+intervals, end = end_followup),by=eval(as.character(x$names$id))]
    grid[,interval:=0:(length(intervals)-1),by=eval(as.character(x$names$id))]
    grid <- pop[,.SD,.SDcols = c(x$names$id)][grid,on = eval(as.character(x$names$id))]
    # FIXME: this does not look great 
    length_interval=unique(round(diff(intervals),0))
    # now awkwardly reset the names
    if (length(x$names$competing)>0 & length(x$long_data$competing_data)>0){
        setnames(x$long_data$competing_data,"competing_date","date")
    }
    if (length(x$names$censoring)>0 & length(x$long_data$censored_data)>0){
        setnames(x$long_data$censored_data,"censored_date","date")
    }
    setnames(x$long_data$outcome_data,"outcome_date","date")
    # FIXME: check for data after the end of followup?
    ## grid <- grid[date<=end+length_interval]
    # mapping outcome information to discrete time scale
    x$data$outcome_data <- widen_outcome(x,
                                         outcome_name = x$names$outcome,
                                         outcome_data = x$long_data$outcome_data,
                                         competing_data = x$long_data$competing_data,
                                         censored_data = x$long_data$censored_data,
                                         grid = grid,
                                         fun_aggregate = NULL,
                                         id = x$names$id)
    # time dependent variables including treatment
    for (Vname in Vnames){
        if (is.list(fun)){
            if (Vname %in% names(fun)){
                # use variable specific function
                vfun <- fun[[Vname]]
            } else{
                # use default
                vfun <- function(x){1*(sum(x)>0)}
            }
        } else {
            vfun <- fun
        }
        stopifnot(is.function(vfun))
        x$data$timevar_data[[Vname]] <- map_grid(grid=grid,
                                                 data=x$long_data$timevar_data[[Vname]],
                                                 name=Vname,
                                                 fun_aggregate = vfun,
                                                 rollforward=(length_interval - 1),
                                                 id = x$names$id)
    }
    return(x)
}



######################################################################
### long_to_wide.R ends here
