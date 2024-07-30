### add_long_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 25 2024 (11:24) 
## Version: 
## Last-Updated: Jul 26 2024 (13:01) 
##           By: Thomas Alexander Gerds
##     Update #: 9
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @export
"add_long_data<-" <- function(x,...,value){
    x$long_data <- copy(value)
    x
}


######################################################################
### add_long_data.R ends here
