### long_to_wide.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 22 2024 (14:07) 
## Version: 
## Last-Updated: Sep 22 2024 (15:30) 
##           By: Thomas Alexander Gerds
##     Update #: 10
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
long_to_wide <- function(x,
                         intervals,
                         fun = function(x){1*(sum(x)>0)}){
    if (length(x$long_data) == 0) {return(NULL)}
    Vnames <- names(x$long_data$timevar_data)
    if (any(duplicated(Vnames))) stop("Duplicated names found in names(x$long_data$timevar_data). Variables must have distinct names.")
    if (missing(intervals)) {stop("Need time intervals to map the long data.")}
    if (length(x$data$baseline_data) == 0 || any(!match(c(x$name_id,"start_followup_date"),names(x$data$baseline_data),nomatch = FALSE)))
        stop("Need baseline data stored as 'x$data$baseline_data' with subject id variable and 'start_followup_date' variable")
    pop <- x$data$baseline_data[,c(x$name_id,"start_followup_date"),with = FALSE]
    #
    # outcome, censored and competing risk define end-of-followup
    #
    if (length(x$long_data$outcome_data)>0)
        if (!("competingrisk_date" %in% names(x$long_data$competingrisk_data)))
            setnames(x$long_data$competingrisk_data,"date","competingrisk_date")
    pop <- x$long_data$competingrisk_data[pop,on = x$name_id]
    if (!("censored_date" %in% names(x$long_data$censored_data)))
        setnames(x$long_data$censored_data,"date","censored_date")
    pop <- x$long_data$censored_data[pop,on = x$name_id]
    if (!("outcome_date" %in% names(x$long_data$outcome_data)))
        setnames(x$long_data$outcome_data,"date","outcome_date")
    pop <- x$long_data$outcome_data[pop,on = x$name_id]
    pop[,end_followup := pmin(censored_date,competingrisk_date,outcome_date,na.rm = TRUE)]
    if (any(is.na(pop$end_followup)))stop("Missing values in end of followup information")
    grid <- pop[,.(date=start_followup_date+intervals, end = end_followup),by=eval(as.character(x$name_id))]
    grid[,interval:=0:(length(intervals)-1),by=eval(as.character(x$name_id))]
    grid <- pop[,.SD,.SDcols = c(x$name_id)][grid,on = eval(as.character(x$name_id))]
    length_interval=unique(round(diff(intervals),0))
    # now awkwardly reset the names
    setnames(x$long_data$competingrisk_data,"competingrisk_date","date")
    setnames(x$long_data$censored_data,"censored_date","date")
    setnames(x$long_data$outcome_data,"outcome_date","date")
    # FIXME: check for data after the end of followup?
    ## grid <- grid[date<=end+length_interval]
    # mapping outcome information to discrete time scale
    x$data$outcome_data <- widen_outcome(outcome_name = x$name_outcome,outcome_data = x$long_data$outcome_data,competingrisk_data = x$long_data$competingrisk_data,censored_data = x$long_data$censored_data,grid = grid,fun.aggregate = NULL,id = x$name_id)
    # time dependent variables including treatment
    for (Vname in Vnames){
        x$data$timevar_data[[Vname]] <- map_grid(grid=grid,data=x$long_data$timevar_data[[Vname]],name=Vname,fun.aggregate = fun,rollforward=(length_interval - 1),id = x$name_id)
    }
    return(x)
}



######################################################################
### long_to_wide.R ends here
