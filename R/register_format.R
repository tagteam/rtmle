### register_format.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: maj  2 2026 (07:14) 
## Version: 
## Last-Updated: sep 20 2026 (06:25) 
##           By: Thomas Alexander Gerds
##     Update #: 17
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
#' @param treatment_variables Optional character vector naming treatment
#'   variables. These variables are returned as exposure intervals with the
#'   columns \code{id}, \code{start_date}, and \code{end_date}; rows where the
#'   treatment is equal to 1 are combined into contiguous intervals. If
#'   omitted, all time-varying variables are returned in the usual
#'   \code{id}, \code{date}, \code{value} format.
#' @return A list with two components:
#'   \itemize{
#'     \item \code{baseline_data}: A data.table containing baseline variables for each individual.
#'     \item \code{timevar_data}: A list of data.tables representing
#'       time-varying information, including event data and visit-related
#'       measurements. Treatment variables named in
#'       \code{treatment_variables} have columns \code{id},
#'       \code{start_date}, and \code{end_date}; other visit variables have
#'       columns \code{id}, \code{date}, and \code{value}.
#'   }
#' @details
#' The function assumes the input object has attributes specifying baseline variables, intermediate events, 
#' absorbing events, baseline visits, visit events, and visit measurements. It extracts relevant information 
#' for each category and organizes it into a structured format. For a
#' treatment variable, a value of 1 denotes exposure. Exposure starts at the
#' corresponding baseline or visit time and ends at the next baseline or visit
#' at which the treatment is no longer 1, or at the subject's last observed
#' time if no later treatment observation is available.
#'
#' @examples
#' data(simulated_cohort)
#' result <- register_format(simulated_cohort,
#'                           treatment_variables = c("A", "B"))
#'
#' @export
register_format <- function(
                            cohort,
                            treatment_variables = NULL
                            ){
    id <- event <- last <- end_date <- start_date <- value <- run <- NULL
    stopifnot(inherits(cohort,"simulated_cohort"))
    info <- attr(cohort,"call",exact = TRUE)
    bvars <- names(eval(info$baseline_variables))
    data.table::setorder(cohort, id, time)
    first <- cohort[, as.integer(.I == .I[1]), by = id][[2]]
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
    if (is.null(treatment_variables)) {
        treatment_variables <- character()
    } else if (!is.character(treatment_variables) ||
               anyNA(treatment_variables) ||
               any(!nzchar(treatment_variables))) {
        stop("Argument treatment_variables must be a character vector of variable names.")
    }
    unknown_treatment_variables <- setdiff(treatment_variables, tvars)
    if (length(unknown_treatment_variables) > 0L) {
        stop("The following treatment variables are not visit variables in the cohort: ",
             paste(unknown_treatment_variables, collapse = ", "))
    }
    treatment_variables <- unique(treatment_variables)
    followup_end <- cohort[, list(end_date = time[.N]), by = id]
    tdata <- lapply(tvars,function(tv){
        if (!(tv %chin% treatment_variables)) {
            out <- cohort[event %chin% c("baseline","visit"),
                          c("id", "time", tv),
                          with = FALSE]
            data.table::setnames(out, c("id", "time", tv),
                                 c("id", "date", "value"))
            return(out)
        }

        treatment_history <- cohort[event %chin% c("baseline","visit"),
                                    c("id", "time", tv),
                                    with = FALSE]
        data.table::setnames(treatment_history,
                             c("id", "time", tv),
                             c("id", "date", "value"))
        data.table::setorder(treatment_history, id, date)
        treatment_history[, run := data.table::rleid(value), by = id]
        treatment_history[, end_date := data.table::shift(date, type = "lead"),
                          by = id]

        ## A treatment value of 1 denotes exposure. Coercing to character
        ## also handles factor-valued treatment columns without relying on
        ## factor/numeric comparison rules.
        exposed <- as.character(treatment_history[["value"]]) %chin% c("1", "TRUE")
        out <- treatment_history[exposed,
                                 list(start_date = first(date),
                                      end_date = last(end_date)),
                                 by = list(id, run)]
        if (nrow(out) > 0L && anyNA(out[["end_date"]])) {
            out[is.na(end_date), end_date :=
                                     followup_end[["end_date"]][match(id, followup_end[["id"]])]]
        }
        out[, run := NULL]
        out[, list(id, start_date, end_date)]
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
