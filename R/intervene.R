### intervene.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: Jul 29 2024 (14:51) 
##           By: Thomas Alexander Gerds
##     Update #: 10
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
intervene <- function(formula,
                      data,
                      protocol,
                      time){
    ## interdata <- copy(model.frame(delete.response(terms(formula(formula))),data,na.action = "na.fail"))
    av <- all.vars(delete.response(terms(formula(formula))))
    if (length(av)>0){
        interdata <- data[,av,with = FALSE]
        ## FIXME: for (tar in targets)
        for (k in 1:nrow(protocol)){
            set(interdata,
                j = protocol[k,variable],
                value = protocol[k,value])
        }
    }else{
        interdata <- NULL
    }
    interdata
}
######################################################################
### intervene.R ends here
