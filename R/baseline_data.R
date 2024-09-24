### baseline_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (08:54) 
## Version: 
## Last-Updated: Sep 22 2024 (15:05) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#'@export
"baseline_data<-" <- function(x,...,value){
    x$data$baseline_data <- value
    x
}



######################################################################
### baseline_data.R ends here
