### long_to_wide.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 22 2024 (14:07) 
## Version: 
## Last-Updated: Apr 10 2025 (11:06) 
##           By: Thomas Alexander Gerds
##     Update #: 37
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
#' This function is used to prepare the discretized time analysis of register data.
#' The start_followup_date is the calendar date where the follow up starts
#' or can be zero if the followup data are readily on the time on study scale  
#' @param x object of class \code{rtmle} 
#' @param intervals a vector of time points that discretize the followup time into intervals
#' @param fun function used to map information onto the intervals. see Details.
#' @details TODO
#' @export
long_to_wide <- function(x,
                         intervals,
                         fun = function(x){1*(sum(x)>0)}){
    start_followup_date = interval = end_followup = censored_date =  competing_date = outcome_date = NULL
    if (length(x$long_data) == 0) {return(NULL)}
    Vnames <- names(x$long_data$timevar_data)
    if (any(duplicated(Vnames))) stop("Duplicated names found in names(x$long_data$timevar_data). Variables must have distinct names.")
    if (missing(intervals)) {stop("Need time intervals to map the long data.")}
    if (length(x$data$baseline_data) == 0 || any(!match(c(x$names$id,"start_followup_date"),names(x$data$baseline_data),nomatch = FALSE)))
        stop("Need baseline data stored as 'x$data$baseline_data' with subject id variable and 'start_followup_date' variable")
    pop <- x$data$baseline_data[,c(x$names$id,"start_followup_date"),with = FALSE]
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
                                         fun.aggregate = NULL,
                                         id = x$names$id)
    # time dependent variables including treatment
    for (Vname in Vnames){
        x$data$timevar_data[[Vname]] <- map_grid(grid=grid,
                                                 data=x$long_data$timevar_data[[Vname]],
                                                 name=Vname,
                                                 fun.aggregate = fun,
                                                 rollforward=(length_interval - 1),
                                                 id = x$names$id)
    }
    return(x)
}



######################################################################
### long_to_wide.R ends here
