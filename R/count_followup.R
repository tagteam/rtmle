### count_followup.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar 24 2025 (15:06)
## Version:
## Last-Updated: Jul 21 2025 (16:03) 
##           By: Thomas Alexander Gerds
##     Update #: 5
#----------------------------------------------------------------------
##
### Commentary:
## This function counts the number of intervals in which the individuals
## are still followed, i.e., at risk of an event. They are followed if they are
## uncensored AND free of outcome AND free of competing events.
##
## The function returns a data.table with person id and the last_interval for
## which the individuals are still followed. E.g.,
##
##  id | last_interval | meaning
##  ----------------------------------------------------------------------
##  1  |      7        | subject 1 is not at risk for intervals 8,9,10, ...
##  2  |      1        | subject 2 is not at risk for intervals 2,3,4, ...
##  3  |      0        | subject 3 is not at risk for intervals 1,2,3, ...
##
## Note that all individuals are
## followed (at-risk) for the first interval.
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
count_followup <- function(x,
                           prepared_data,
                           outcome_variables,
                           censoring_variables,
                           competing_variables,
                           max_time_horizon){
    x$followup <- prepared_data[,c(x$names$id),with = FALSE]
    set(x$followup,j = "last_interval",value = numeric(NROW(x$followup)))
    if (max_time_horizon>1){
        for (j in 0:(max_time_horizon-2)){
              vital <- 1*(prepared_data[[outcome_variables[[j+1]]]] %in% "0")
            if (length(censoring_variables)>0){
                vital <- vital * 1*(prepared_data[[censoring_variables[[j+1]]]] %in% x$names$uncensored_label)
            }
            if (length(competing_variables)>0)
                vital <- vital * (prepared_data[[competing_variables[[j+1]]]] %in% "0")
            set(x$followup,j = "last_interval",value = x$followup[["last_interval"]]+vital)
        }
    }
    x
}


######################################################################
### count_followup.R ends here
