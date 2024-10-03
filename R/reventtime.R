### reventtime.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 12 2024 (09:38) 
## Version: 
## Last-Updated: Oct  2 2024 (16:03) 
##           By: Thomas Alexander Gerds
##     Update #: 69
#----------------------------------------------------------------------
## 
### Comaxtimeentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
reventtime <- function(n,
                       breaks,
                       cumhazard,
                       hazardratio,
                       entrytime = NULL,
                       decimals = NULL){
    if (missing(n)) stop("Please specify sample size 'n'")
    if (any(hazardratio == 0)) stop("Hazardratios cannot contain zeros")
    if (all(breaks>0)) {
        breaks <- c(0,breaks)
        cumhazard <- c(0,cumhazard)
    }
    ## linear approximation
    lapprox <- function(values, breaks, fx){
        pos <- prodlim::sindex(jump.times = breaks, eval.times = values)
        maxindex <- which(pos == length(breaks))
        next_pos <- pos + 1
        pos <- pmax(pos,1)
        next_pos[maxindex] <- length(breaks)
        approx_value <- (values - breaks[pos])/(breaks[next_pos] - breaks[pos])
        approx_value[maxindex] <- 0
        res <- approx_value * (fx[next_pos] - fx[pos]) + fx[pos]
        res[is.na(res)] <- utils::tail(fx, 1)
        return(res)
    }
    if (missing(hazardratio)) {
        hazardratio <- rep(1,n)
    }
    maxtime <- utils::tail(breaks, 1)
    if (!is.null(entrytime) && any(entrytime>0)) {
        if (length(entrytime) == 1) entrytime <- rep(entrytime, n)
        ## if (FALSE){
        pos <- prodlim::sindex(jump.times = breaks, eval.times = entrytime)
        maxindex <- which(pos == length(breaks))
        next_pos <- pos + 1
        pos <- pmax(pos,1)
        next_pos[maxindex] <- length(breaks)
        approx_value <- (entrytime - breaks[pos])/(breaks[next_pos] - breaks[pos])
        approx_value[maxindex] <- 0
        entry_cumhazard <- approx_value * (cumhazard[next_pos] - cumhazard[pos]) + cumhazard[pos]
        entry_cumhazard[is.na(entry_cumhazard)] <- utils::tail(cumhazard, 1)
        ## }
        ## entry_cumhazard <- lapprox(values = entrytime,breaks = breaks,fx = cumhazard)
    }
    else {
        entry_cumhazard <- rep(0, n)
    }
    # The distribution of T|X is the same as Lambda^-1(E/HR)
    # where E is expontential with rate 1
    erate <- stats::rexp(n)/hazardratio + entry_cumhazard
    ## if (FALSE){
    pos <- prodlim::sindex(jump.times = cumhazard, eval.times = erate)
    maxindex <- which(pos == length(cumhazard))
    next_pos <- pos + 1
    pos <- pmax(pos,1)
    next_pos[maxindex] <- length(cumhazard)
    approx_value <- (erate - cumhazard[pos])/(cumhazard[next_pos] - cumhazard[pos])
    approx_value[maxindex] <- 0
    etime <- approx_value * (breaks[next_pos] - breaks[pos]) + breaks[pos]
    etime[is.na(etime)] <- utils::tail(breaks, 1)
    ## }
    ## etime <- pmin(lapprox(values = erate,fx = breaks,breaks = cumhazard),maxtime)
    if (!is.null(decimals))
        etime = round(etime,decimals)
    return(etime)
}


######################################################################
### reventtime.R ends here
