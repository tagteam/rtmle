### intervene.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: Oct 16 2024 (15:37) 
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
intervene <- function(data,
                      intervention_table,
                      time){
    ## FIXME: is it optional to set variables at all times or only the current and later?
    ##        believe not but need to check Schnitzer et al.
    interdata <- copy(data)
    for (k in 1:nrow(intervention_table)){
        set(interdata,
            j = intervention_table[k][["variable"]],
            value = intervention_table[k][["value"]])
    }
    interdata
}
######################################################################
### intervene.R ends here
