widen_outcome <- function(x,
                          outcome_data = NULL,
                          censored_data = NULL,
                          competing_data = NULL,
                          grid = NULL,
                          fun_aggregate = NULL){
    # The explicit event tables are prepared by long_to_wide(). Allow the
    # outcome table to be omitted for direct internal calls; the normal path
    # passes all event tables after applying the follow-up precedence rules.
    if (is.null(outcome_data)) {
        if (length(x$long_data$outcome_data) == 0) {
            stop("widen_outcome: Object does not contain outcome data at x$long_data$outcome_data")
        }
        outcome_data <- data.table::copy(x$long_data$outcome_data)
        if ("date" %in% names(outcome_data)) {
            data.table::setnames(outcome_data, "date", "outcome_date")
        }
    }
    if (is.null(grid)) stop("widen_outcome: a grid is required.")
    # -----------------------------------------------------------------------
    # death and right censored
    # -----------------------------------------------------------------------
    #
    # Notes:
    #       a) when outcome or death has occurred the value 1 persists, i.e.,
    #          the last observation is carried forward.
    #          this is done by discretize.
    #       b) when outcome occurs before death or censored then
    #          the value of death or censored was removed before calling discretize.
    #       c) once censored both outcome and death variables are NA
    # 
    censored_variables <- NULL
    if (length(x$names$censoring)>0 && !is.null(censored_data) &&
        NROW(censored_data)>0){
        if (any(duplicated(censored_data[[x$names$id]]))){
            stop("Duplicated person id in dates for censoring risks.")
        }
        current_censored_data <- data.table::copy(censored_data)
        if ("censored_date" %in% names(current_censored_data)) {
            data.table::setnames(current_censored_data, "censored_date", "date")
        }
        censored_variables <- discretize(
            method = "event",
            data = current_censored_data,
            grid = grid,
            name = x$names$censoring,
            id = x$names$id,
            threshold = NULL,
            lookback_window = Inf,
            # the order must be censored_label, uncensored_label
            values = rev(x$names$censored_levels),
            fun_aggregate = fun_aggregate,
            fill = NA
        )
        # This makes sure that all censored variables are factors with levels
        # ordered as c(uncensored, censored).
        for (cc in names(censored_variables)[-1]){
            set(censored_variables,
                j = cc,
                value = factor(censored_variables[[cc]],
                               levels = x$names$censored_levels))
        }
    }
    competing_variables <- NULL
    if (length(x$names$competing)>0 && !is.null(competing_data) &&
        NROW(competing_data)>0){
        if (any(duplicated(competing_data[[x$names$id]]))){
            stop("Duplicated person id in dates for competing risks.")
        }
        current_competing_data <- data.table::copy(competing_data)
        if ("competing_date" %in% names(current_competing_data)) {
            data.table::setnames(current_competing_data, "competing_date", "date")
        }
        competing_variables <- discretize(
            method = "event",
            data = current_competing_data,
            grid = grid,
            name = x$names$competing,
            id = x$names$id,
            threshold = NULL,
            lookback_window = Inf,
            values = c(1, 0),
            fun_aggregate = fun_aggregate,
            fill = NA
        )
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
    current_outcome_data <- data.table::copy(outcome_data)
    if ("outcome_date" %in% names(current_outcome_data)) {
        data.table::setnames(current_outcome_data, "outcome_date", "date")
    }
    outcome_variables <- discretize(
        method = "event",
        data=current_outcome_data,
        grid = grid,
        name=x$names$outcome,
        id = x$names$id,
        threshold = NULL,
        lookback_window = Inf,
        values = c(1,0),
        fun_aggregate = fun_aggregate,
        fill = NA
    )
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
