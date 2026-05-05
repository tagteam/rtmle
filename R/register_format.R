### register_format.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: maj  2 2026 (07:14) 
## Version: 
## Last-Updated: maj  4 2026 (12:02) 
##           By: Thomas Alexander Gerds
##     Update #: 13
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Register Format Conversion for Simulated Cohort Data
#'
#' This function transforms a simulated cohort object into a structured format 
#' containing baseline data, time-varying data, and event data. The output is 
#' a list with separate components for baseline and time-varying information, 
#' facilitating further analysis.
#'
#' @param cohort A data object of class \code{"simulated_cohort"}.
#' @return A list with two components:
#'   \itemize{
#'     \item \code{baseline_data}: A data.table containing baseline variables for each individual.
#'     \item \code{timevar_data}: A list of data.tables representing time-varying information, including event data and visit-related measurements.
#'   }
#' @details
#' The function assumes the input object has attributes specifying baseline variables, intermediate events, 
#' absorbing events, baseline visits, visit events, and visit measurements. It extracts relevant information 
#' for each category and organizes it into a structured format.
#'
#' @examples
#' data(simulated_cohort)
#' result <- register_format(simulated_cohort)
#'
#' @export
register_format <- function(cohort){
    id <- event <- NULL
    stopifnot(inherits(cohort,"simulated_cohort"))
    info <- attr(cohort,"call",exact = TRUE)
    bvars <- names(eval(info$baseline_variables))
    data.table::setorder(cohort, id, time)
    first <- cohort[, as.integer(.I == .I[1]), by = id][[2]]
    last <- cohort[, as.integer(.I == .I[.N]), by = id][[2]]
    baseline_data <- cohort[first == 1L, c("id", bvars),with = FALSE]
    evars <- c(names(eval(info$intermediate_events)),
               names(eval(info$absorbing_events)))
    edata <- lapply(evars,function(e){
        cohort[event == e, list(id, date = time)]
    })
    names(edata) <- evars
    tvars <- unique(c(names(eval(info$baseline_visit)),
                      names(eval(info$visit_events)),
                      names(eval(info$visit_measurements))))
    tdata <- lapply(tvars,function(tv){
        setnames(cohort[event %chin% c("baseline","visit"), c("id", "time", tv),with = FALSE],c("id","date","value"))
    })
    names(tdata) <- tvars
    reg_format <- list(
        baseline_data = baseline_data[],
        timevar_data = c(tdata,edata)
    )
    reg_format
}

######################################################################
### register_format.R ends here
