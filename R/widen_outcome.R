widen_outcome <- function(x,
                          outcome_name,
                          outcome_data,
                          competing_data,
                          censored_data,
                          grid,
                          fun_aggregate = NULL,
                          id){
    stopifnot(id %in% names(outcome_data))
    # -----------------------------------------------------------------------
    # death and right censored
    # -----------------------------------------------------------------------
    if (length(x$names$competing)>0 & length(x$long_data$competing_data)>0){
        wide <- map_grid(grid=grid,
                         data=competing_data,
                         name="Dead",
                         rollforward=Inf,
                         id = id)
    }else{
        wide <- NULL
    }
    if (length(x$names$censoring)>0 & length(x$long_data$censored_data)>0){
        # naturally, NA values mean censored. hence: fill=1
        w <- map_grid(grid=grid,
                   data=censored_data,
                   name="Censored",
                   rollforward=Inf,
                   values=c("censored","uncensored"),
                   fun_aggregate = fun_aggregate,
                   fill="censored",
                   value_is_factor=TRUE,
                   id = id)
        if (length(wide)>0) wide <- wide[w] else wide <- w
    }
    # -----------------------------------------------------------------------
    # only interested in new outcomes with onset after index
    # but want to tag patients who are in hospital with the outcome
    # at the index date, in order to use this as a baseline variable
    ## hospital diagnoses overlapping start
    ## admitted_index=outcome_data[date<=start & discharge>start,unique(id)]
    ## outcome_data=outcome_data[date>start]
    ## only interested in first new outcome
    ## outcome_data=outcome_data[outcome_data[,.I[1],by=id]$V1]
    w <- map_grid(grid=grid,
                  data=outcome_data,
                  name=outcome_name,
                  rollforward=Inf,
                  id = id)
    # when outcome occurs at baseline re-set to value 2
    ## set(w,i=which(wide$id%in%admitted_index),j=paste0(outcome,"_0"),value=2)
    if (length(wide)>0) wide <- wide[w] else wide <- w
    #
    # Notes:
    #       a) when outcome or death has occurred the value 1 persists, i.e.,
    #          the last observation is carried forward
    #       b) outcome in the interval which starts at index does not count as outcome!
    #       c) patients who die at index date are excluded
    #          hence the loop starts at index 1
    #       d) when outcome occurs before death/censored in same interval
    #          then the value of death/censored does not matter
    #          exception: outcome, death, censored at index date
    #       e) when censoring occurs but not outcome then outcome is NA
    max_interval = max(grid$interval)
    for (tk in 1:(max_interval-1)){
        Ok=paste0(outcome_name,"_",tk)
        Ok_next=paste0(outcome_name,"_",tk+1)
        Dk=paste0("Dead","_",tk)
        Dk_next=paste0("Dead","_",tk+1)
        Ck=paste0("Censored","_",tk)
        Ck_next=paste0("Censored","_",tk+1)
        # when died at t_k then also at t_{k+1}
        if (any(has_died <- wide[[Dk]]==1))
            set(wide,j=Dk_next,i=which(has_died),value=1)
        # when censored at t_k then also at t_{k+1}
        #   --- nothing to do ---
        # when outcome at t_k then also at t_{k+1}
        if (any(has_outcome <- wide[[Ok]]==1)){
            set(wide,j=Ok_next,i=which(has_outcome),value=1)
        }
        # last outcome value carried forward
        if (any(miss_outcome <- is.na(wide[[Ok_next]])))
            set(wide,j=Ok_next,i=which(miss_outcome),value=wide[miss_outcome][[Ok]])
        # when censored in interval but not outcome, then outcome and death are both NA
        has_censored <- wide[[Ck]]==1
        if (any(has_censored & !(has_outcome)))
            set(wide,j=Ok_next,i=which(has_censored & !(has_outcome)),value=NA)
        if (any(has_censored & !(has_died)))
            set(wide,j=Dk_next,i=which(has_censored & !(has_died)),value=NA)
    }
    wide
}
