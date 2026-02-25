widen_outcome <- function(x,
                          grid,
                          fun_aggregate = NULL){
    if (length(x$long_data$outcome_data) == 0) {
        stop("widen_outcome: Object does not contain outcome data at x$long_data$outcome_data")
    }
    # -----------------------------------------------------------------------
    # death and right censored
    # -----------------------------------------------------------------------
    #
    # Notes:
    #       a) when outcome or death has occurred the value 1 persists, i.e.,
    #          the last observation is carried forward.
    #          this is done by map_grid.
    #       b) when outcome occurs before death or censored then
    #          the value of death or censored was removed before calling map_grid
    #       c) once censored both outcome and death variables are NA
    # 
    censored_variables <- NULL
    if (length(x$names$censoring)>0 && length(x$long_data$censored_data)>0){
        if (any(duplicated(x$long_data$censored_data[[x$names$id]]))){
            stop("Duplicated person id in dates for censoring risks.")
        }else{
            censored_before_outcome_and_death <- data.table(id = setdiff(x$long_data$censored_data[[x$names$id]],
                                                                         c(x$long_data$outcome_data[[x$names$id]],
                                                                           x$long_data$competing_data[[x$names$id]])))
            data.table::setnames(censored_before_outcome_and_death,x$names$id)
            data.table::setkeyv(censored_before_outcome_and_death,x$names$id)
        }
        if (NROW(censored_before_outcome_and_death)>0){
            censored_variables <- map_grid(grid=grid,
                                           data=x$long_data$censored_data[censored_before_outcome_and_death,on = x$names$id],
                                           name=x$names$censoring,
                                           rollforward=Inf,
                                           # the order must be censored_label, uncensored_label 
                                           values=rev(x$names$censored_levels),
                                           fun_aggregate = fun_aggregate,
                                           id = x$names$id)
            # this makes sure that all censored variables are
            # factors with levels order as c(uncensored,censored)
            for (cc in names(censored_variables)[-1]){
                set(censored_variables,j=cc,value=factor(censored_variables[[cc]],levels=x$names$censored_levels))
            }
        }
    }
    competing_variables <- NULL
    if (length(x$names$competing)>0 && length(x$long_data$competing_data)>0){
        if (any(duplicated(x$long_data$competing_data[[x$names$id]]))){
            stop("Duplicated person id in dates for competing risks.")
        }else{
            # when both outcome and competing risk occurs then the competing risk date is ignored 
            competing_risk_before_outcome <- data.table(id = setdiff(x$long_data$competing_data[[x$names$id]],
                                                                     x$long_data$outcome_data[[x$names$id]]))
            data.table::setnames(competing_risk_before_outcome,x$names$id)
            data.table::setkeyv(competing_risk_before_outcome,x$names$id)
        }
        if (NROW(competing_risk_before_outcome)>0){
            competing_variables <- map_grid(grid=grid,
                                            data=x$long_data$competing_data[competing_risk_before_outcome,on = x$names$id],
                                            name=x$name$competing,
                                            rollforward=Inf,
                                            id = x$names$id)
        }
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
    outcome_variables <- map_grid(grid=grid,
                                  data=x$long_data$outcome_data,
                                  name=x$names$outcome,
                                  rollforward=Inf,
                                  id = x$names$id)
    # Once censoring has occurred all following competing and outcome variables
    # should be NA. Note that by construction id's where both censored dates AND
    # outcome/competing dates are available the censored dates are removed before
    # at the call of map_grid
    if (length(censored_variables)>0){
        cens_varnames <- setdiff(names(censored_variables),x$names$id)
        # outcome
        out_varnames <- setdiff(names(outcome_variables),x$names$id)
        has_competing <- length(competing_variables)>0
        if (has_competing){
            comp_varnames <- setdiff(names(competing_variables),x$names$id)
        }
        for (j in 1:length(cens_varnames)){
            # allow missing values in censored_variables
            has_censored <- censored_variables[[cens_varnames[[j]]]] %in% x$names$censored_label
            if (any(has_censored)){
                set(outcome_variables,j = out_varnames[[j]], i = which(has_censored),value = NA)
                if (has_competing){
                    set(competing_variables,j = comp_varnames[[j]], i = which(has_censored),value = NA)
                }
            }
        }
        # join on id
        wide <- outcome_variables[censored_variables]
        # competing
        if (has_competing){
            # join on id
            wide <- wide[competing_variables]
        }
    }else{
        if (length(competing_variables)>0){
            # join on id
            wide <- competing_variables[outcome_variables]
        }else{
            wide <- outcome_variables
        }
    }
    wide
}
