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
#' This is the default function used by \code{\link{protocol}}. It is also a
#' convenient starting point inside a user-defined intervention function:
#' first apply the nominal protocol values, then replace treatment values for
#' rows governed by a history-dependent rule.
#'
#' @param data The current history data.
#' @param intervention_table Long-format intervention table supplied by
#'   \code{protocol()}.
#' @param time_node Current intervention node. Only treatment assignments up to
#'   and including this node are applied.
#' @return A copy of \code{data} with protocol treatment values applied.
#' @seealso \code{\link{protocol}}
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
        set(interdata,
            j = intervention_table[k][["variable"]],
            value = rep(intervention_table[k][["value"]],N))
    }
    interdata
}
######################################################################
### intervene.R ends here
