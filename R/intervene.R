### intervene.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: Mar  7 2025 (18:41) 
##           By: Thomas Alexander Gerds
##     Update #: 26
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
    N <- NROW(interdata)
    for (k in 1:nrow(intervention_table)){
        set(interdata,
            j = intervention_table[k][["variable"]],
            value = rep(intervention_table[k][["value"]],N))
    }
    interdata
}
######################################################################
### intervene.R ends here
