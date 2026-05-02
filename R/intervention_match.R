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
##' Check adherence to a protocol
##'
##' Builds the \code{intervention_match} matrix for a protocol. The function is
##' called by \code{\link{protocol}} when \code{\link{prepare_rtmle_data}} has
##' already been run, and by \code{\link{run_rtmle}} if the matrix is still
##' missing at estimation time.
##'
##' @title Check adherence to a protocol
##' @param x An \code{rtmle} object as returned by \code{\link{rtmle_init}}.
##' @param protocol_name Name of the protocol to check.
##' @return The modified \code{rtmle} object.
##' @seealso \code{\link{protocol}}, \code{\link{prepare_rtmle_data}},
##'   \code{\link{run_rtmle}}, \code{\link{plot_adherence}}
##' @examples
##' x <- rtmle_init(time_grid = 0:2, name_id = "id", name_outcome = "Y")
##' x$prepared_data <- data.table::data.table(
##'     id = 1:4,
##'     A_0 = factor(c("1", "1", "0", "1"), levels = c("0", "1")),
##'     A_1 = factor(c("1", "0", "0", "1"), levels = c("0", "1")))
##' x <- protocol(x, name = "Always_A",
##'               intervention = data.frame(time = x$intervention_nodes,
##'                                         A = factor("1", levels = c("0", "1"))))
##' x$protocols$Always_A$intervention_match <- NULL
##' x <- intervention_match(x, "Always_A")
##' x$protocols$Always_A$intervention_match
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
