### treatment_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 16 2024 (10:28) 
## Version: 
## Last-Updated: Sep 16 2024 (10:29) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

"treatment_data<-" <- function(x,...,value){
    x$long_data$treatment_data <- value
    x
}

######################################################################
### treatment_data.R ends here
