### intervene.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: Nov  1 2024 (10:47) 
##           By: Thomas Alexander Gerds
##     Update #: 21
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
intervene <- function(data,
                      intervention_table,
                      time){
    ## FIXME: is it optional to set variables at all times or only the current and later?
    ##        believe not but need to check Schnitzer et al.
    interdata <- copy(data)
    for (k in 1:nrow(intervention_table)){
        if (is.character(intervention_table[k][["value"]]) && is.function(eval(as.name(intervention_table[k][["value"]])))){
            set(interdata,
                j = intervention_table[k][["variable"]],
                value = do.call(intervention_table[k][["value"]],list(data = interdata,time = k)))
        } else{
            set(interdata,
                j = intervention_table[k][["variable"]],
                value = intervention_table[k][["value"]])
        }
    }
    interdata
}
######################################################################
### intervene.R ends here
