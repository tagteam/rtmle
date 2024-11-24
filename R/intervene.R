### intervene.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: Nov 24 2024 (06:50) 
##           By: Thomas Alexander Gerds
##     Update #: 23
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
