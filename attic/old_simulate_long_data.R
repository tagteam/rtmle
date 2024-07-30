### old_simulate_long_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 11 2024 (13:24) 
## Version: 
## Last-Updated: Jul 17 2024 (10:49) 
##           By: Thomas Alexander Gerds
##     Update #: 55
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#'
#' library(data.table)
#' library(mets)
#' library(ggplot2)
#' # baseline hazard functions 
#' d <- data.table(breaks = seq(1,10*365,1))
#' d[,cumhaz_A := cumsum(rep(0.002,.N))]
#' d[,cumhaz_L := cumsum(rep(0.001,.N))]
#' d[,cumhaz_Y := cumsum(rep(0.001,.N))]
#' d[,cumhaz_D := cumsum(rep(0.001,.N))]
#' d[,cumhaz_C := cumsum(rep(0.001,.N))]
#' ## ggplot(d,aes(x = breaks,y = cumhaz_Y))+geom_line()+geom_line(aes(x = breaks,y = cumhaz_D,color = 2))
#' set.seed(18)
#' n <- 20
#' rf <- list(A = rep(20,n),L = rep(10,n),Y = rep(1,n),D = rep(1,n),C = rep(1,n))
#' ld <- old_simulate_long_data(n,breaks = d$breaks,baseline_hazard = list(A = d$cumhaz_A,L = d$cumhaz_L,Y = d$cumhaz_Y,D = d$cumhaz_D,C = d$cumhaz_C),hazardratio = rf)
#' ld[,table(table(id))]
#' ld
old_simulate_long_data <- function(n,
                               breaks,
                               baseline_hazard,
                               hazardratio,
                               max.recurrent = 100,
                               var.z = 0.22,
                               cor.mat = NULL, 
                               ...) {
    if (missing(n)){
        if (missing(hazardratio)) stop("Specify Either 'n' or 'hazardratio'")
        n <- length(hazardratio$A)
    } else{
        if (missing(hazardratio) || length(hazardratio) == 0)
            hazardratio <- list(A = rep(1,n),L = rep(1,n),Y = rep(1,n),D = rep(1,n),C = rep(1,n))
    }
    status <- survival_status <- survival_time <- NULL
    # frailty
    zz <- rgamma(n, 1/var.z[1]) * var.z[1]
    # terminal events and right censored
    tall <- data.table(do.call(cbind,lapply(c("D","Y","C"),function(name){
        x = reventtime(breaks = breaks,cumhazard = baseline_hazard[[name]],hazardratio = (zz * hazardratio[[name]]),entrytime = 0,decimals = 2)
        x[]
    })))
    setnames(tall,paste0("time_",c("D","Y","C")))
    tall[,survival_time := do.call(pmin,.SD),.SDcols = intersect(c("time_D","time_C"),names(tall))]
    tall[,Y_time := do.call(pmin,.SD),.SDcols = intersect(c("time_Y","time_D","time_C"),names(tall))]
    tall[,survival_status := 0]
    tall[time_D <= time_C,survival_status := 1]
    for (name in intersect(names(tall),c("time_D","time_Y","time_C")))
        set(tall,j = name,value = NULL)
    tall[,Y_status := 2*survival_status]
    tall[Y_time < survival_time,Y_status := survival_status]
    tall[, id := 1:n]
    still <- copy(tall)
    still[,event_time := 0]
    i <- 1
    while (i < max.recurrent && nrow(still)>0){
        i <- i + 1
        # recurrent L events
        still <- still[,data.table(id = id,survival_time,Y_time,survival_status,Y_status,event_type = "A",event_time = reventtime(breaks = breaks,cumhazard = baseline_hazard[["A"]],hazardratio = (zz[id] * hazardratio[["A"]][id]),entrytime = event_time,decimals = 2))[Y_time>event_time]]
        tall <- rbindlist(list(tall,still),fill = TRUE)
        if (nrow(still)>0){
            still <- still[,{
                data.table(id = id,survival_time,Y_time,survival_status,Y_status,event_type = "L",event_time = reventtime(breaks = breaks,cumhazard = baseline_hazard[["L"]],hazardratio = (zz[id] * hazardratio[["L"]][id]),entrytime = event_time,decimals = 2))[Y_time>event_time]
            }]
            tall <- rbindlist(list(tall,still),fill = TRUE)
        }
    }
    tall[is.na(event_type),`:=`(event_type="Y",event_time = Y_time)]
    setcolorder(tall,"id")
    setkey(tall,id,event_time)
    tall[]
}



######################################################################
### old_simulate_long_data.R ends here
