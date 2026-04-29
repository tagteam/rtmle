### intervention_match.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: apr 23 2026 (16:09) 
## Version: 
## Last-Updated: apr 29 2026 (07:36) 
##           By: Thomas Alexander Gerds
##     Update #: 13
#----------------------------------------------------------------------
## 
### Commentary: 
#
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Create a matrix called intervention_match 
##'
##' This function is run by \code{protocol} only if \code{prepare_rtmle_data} is run before adding a protocol.
##' It is also run by \code{run_rtmle} if needed at that stage.
##' @title Checking adherence to a given protocol (treatment regimen, intervention).
##' @param x An rtmle object as obtained with \code{rtmle_init}.
##' @param protocol_name Name of the protocol to check.
##' @return The modified object
##' @seealso protocol 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
intervention_match <- function(x,protocol_name){
    variable <- NULL
    N <- NROW(x$prepared_data)
    # define a matrix which indicates if the intervention is followed
    # the matrix should have a row for each individual and a column for
    # each intervention node (time interval)
    #
    intervention_table <- na.omit(x$protocols[[protocol_name]]$intervention_table)
    if (
        N>0 && 
        # no need to redo this operation ever
        length(intervention_match <- x$protocols[[protocol_name]]$intervention_match) == 0
    ){
        intervention_match <- matrix(0,ncol = length(unique(intervention_table$time)),nrow = N)
        intervention_match_names <- vector("character",NCOL(intervention_match))
        previous <- rep(1,N)
        for(k in unique(intervention_table$time)){
            intervention_variables <- intervention_table[time == k][["variable"]]
            if (length(intervention_variables)>0){
                observed_values <- x$prepared_data[,intervention_variables,with = FALSE]
                for (v in 1:length(intervention_variables)){
                    # when there are multiple intervention variables
                    # all observed values must match
                    intervention_values <- intervention_table[time == k & variable == intervention_variables[[v]]][["value"]]
                    intervention_match[,k+1] <- previous <- previous*(observed_values[[intervention_variables[[v]]]] %in% intervention_values)
                    # when there are multiple treatment variables we paste-collapse the names
                }
                intervention_match_names[k+1] <- paste0(intervention_variables,collapse = ",")
            }
            # else{
                # NO INTERVENTION AT THAT NODE
            #}
        }
        ## the intervention_match matrix has one column per time point
        colnames(intervention_match) <- intervention_match_names
        x$protocols[[protocol_name]]$intervention_match <- intervention_match
    }
    x
}


######################################################################
### intervention_match.R ends here
