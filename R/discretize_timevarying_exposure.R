### discretize_exposure.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 24 2026 (11:04) 
## Version: 
## Last-Updated: mar 23 2026 (10:54) 
##           By: Thomas Alexander Gerds
##     Update #: 45
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# data contains exposure episodes
# grid contains intervals on a discrete time scale
# the function discretizes how many days a subject
# was exposed within the interval and defines exposure when the
# number of days exceeds threshold. E.g., when the length of the
# interval is 6 months (180 days) and the threshold is set to 91 then exposure
# is 1 if the treatment duration exceeds 50% of the length interval
# and 0 otherwise
discretize_timevarying_exposure <- function(data,
                                            grid,
                                            name,
                                            point_exposure,
                                            threshold,
                                            id){
    exposure = start_exposure = end_exposure = end_interval = start_interval = interval = NULL
    data <- copy(data)
    grid <- copy(grid)
    if (point_exposure){
        setnames(data,"date","start_exposure")
        set(data,j = "end_exposure",value = data[["start_exposure"]])
    }
    setkeyv(data,c(id, "start_exposure", "end_exposure"))
    setkeyv(grid,c(id, "start_interval", "end_interval"))
    overlap <- data.table::foverlaps(
                               x = grid,
                               by.x = c(id, "start_interval", "end_interval"),
                               y = data,
                               by.y = c(id, "start_exposure", "end_exposure"),
                               type = "any",
                               nomatch = NA
                           )
    if (point_exposure){
        overlap[,exposure :=
                     1*(
                         # baseline exposure 
                         (interval == 0 & start_exposure == 0) |
                         (start_exposure <= end_interval & start_exposure >= start_interval)
                     )]
        # missing values occur in intervals without any overlap
        overlap[is.na(exposure),exposure := 0]
        # summarize across interval
        overlap <- overlap[,list(exposure = 1*(any(exposure == 1))),by = c(id,"interval")]
    }else{
        overlap[,exposure := (pmin(end_interval,end_exposure)-pmax(start_interval,start_exposure))]
        overlap[is.na(exposure),exposure := 0]
        # baseline exposure
        overlap[interval == 0 & start_exposure == 0, exposure := 1]
        # total duration in interval
        setkeyv(overlap, c(id, "interval"))
        overlap <- overlap[, list(exposure = sum(exposure)), by = c(id, "interval")]
        if (is.numeric(threshold)){
            set(overlap,j = "exposure",value = 1*(overlap[["exposure"]]>threshold))
        }
    }
    xx <- fast_cast(overlap, id = id, value_col = "exposure", fill = NA)
    setkey(xx,id)
    ## xx1 <- data.table::dcast(overlap,id~interval,value.var="exposure",sep="_",fill=NA)
    setnames(xx,c(id,paste0(name,"_",names(xx)[-1])))
    xx
}

######################################################################
### discretize_exposure.R ends here
