### intervene.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: maj 21 2026 (08:35) 
##           By: Thomas Alexander Gerds
##     Update #: 27
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Apply the static treatment values in an intervention table
#'
#' This is the default function used by \code{\link{regime}}. It is also a
#' convenient starting point inside a user-defined intervention function:
#' first apply the nominal regime values, then replace treatment values for
#' rows governed by a history-dependent rule.
#'
#' @param data The current history data.
#' @param intervention_table Long-format intervention table supplied by
#'   \code{regime()}.
#' @param time_node Current intervention node. Only treatment assignments up to
#'   and including this node are applied.
#' @return A copy of \code{data} with regime treatment values applied.
#' Binary factor values in the intervention table are coerced to the type of
#' the corresponding data column, so integer treatment histories are handled
#' in the same way as factor histories.
#' @seealso \code{\link{regime}}
#' @export
intervene <- function(data,
                      intervention_table,
                      time_node){
    interdata <- copy(data)
    N <- NROW(interdata)
    if ("time_node" %in% names(intervention_table)) {
        intervention_table <- intervention_table[
            !is.na(intervention_table[["time_node"]]) &
                intervention_table[["time_node"]] <= time_node
        ]
    }
    if ("value" %in% names(intervention_table)) {
        intervention_table <- intervention_table[!is.na(intervention_table[["value"]])]
    }
    for (k in seq_len(nrow(intervention_table))){
        variable <- intervention_table[k][["variable"]]
        value <- intervention_table[k][["value"]]
        # Intervention tables are normally built from factors, whereas data
        # imported from register-format event histories often store binary
        # treatments as integers.  Assigning a factor directly to an integer
        # column would use its internal level code (so factor("0") becomes
        # 1).  Coerce the table value to the representation used by the
        # supplied data before assigning it.
        current <- interdata[[variable]]
        if (is.factor(current)) {
            value <- factor(as.character(value), levels = levels(current))
        } else if (is.factor(value)) {
            value <- as.character(value)
            if (is.integer(current)) {
                value <- as.integer(value)
            } else if (is.numeric(current)) {
                value <- as.numeric(value)
            } else if (is.logical(current)) {
                value <- value %in% c("TRUE", "true", "1")
            }
        }
        set(interdata,
            j = variable,
            value = rep(value,N))
    }
    interdata
}
######################################################################
### intervene.R ends here
