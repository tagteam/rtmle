### intervention_match.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: apr 23 2026 (16:09) 
## Version: 
## Last-Updated: apr 23 2026 (17:44) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
#
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
intervention_match <- function(x,protocol_name){
    variable <- NULL
    N <- NROW(x$prepared_data)
    # define a matrix which indicates if the intervention is followed
    # the matrix should have a row for each individual and a column for
    # each intervention node (time interval)
    #
    if (
        N>0 && 
        # no need to redo this operation ever
        length(intervention_match <- x$protocols[[protocol_name]]$intervention_match) == 0
    ){
        intervention_match <- matrix(0,ncol = length(x$intervention_nodes),nrow = N)
        intervention_match_names <- vector("character",length(x$intervention_nodes))
        previous <- rep(1,N)
        for(k in x$intervention_nodes){
            intervention_variables <- x$protocols[[protocol_name]]$intervention_table[time == k][["variable"]]
            if (length(intervention_variables)>0){
                observed_values <- x$prepared_data[,intervention_variables,with = FALSE]
                for (v in 1:length(intervention_variables)){
                    # when there are multiple intervention variables
                    # all observed values must match
                    intervention_values <- x$protocols[[protocol_name]]$intervention_table[time == k & variable == intervention_variables[[v]]][["value"]]
                    intervention_match[,k+1] <- previous <- previous*(observed_values[[intervention_variables[[v]]]] %in% intervention_values)
                    # when there are multiple treatment variables we paste-collapse the names
                    intervention_match_names[k+1] <- paste0(intervention_variables,collapse = ",")
                }
            }else{
                # FIXME: does this work when we only intervene on baseline treatment but not on subsequent treatments?
                # last value carried forward
                intervention_match[,k+1] <- previous
                intervention_match_names[k+1] <- paste0("No_intervention",k)
            }
        }
        ## the intervention_match matrix has one column per time point
        colnames(intervention_match) <- intervention_match_names
        x$protocols[[protocol_name]]$intervention_match <- intervention_match
    }
    x
}


######################################################################
### intervention_match.R ends here
