## Define a collection of history functions that maps long
## time-varying data to wide format using grid.
##
## TODO: implement more functions
history_method <- function(data, id, grid, name, method){

    ## Just apply method if it is a function
    if(is.function(method)){
        out = method(data = data,
                     grid = grid,
                     name = name,
                     id = id)        
        return(out)
    }

    ## Methods names check:
    method_reqs = list(last = c("date","value"),
                       exposed = c("start_exposure","end_exposure"),
                       exposure_time = c("start_exposure","end_exposure"),
                       event = c("date"))

    ## Attempt to guess method if not specified:
    if(is.null(method)){
        matches <- vapply(method_reqs,
                          function(req) all(req %in% names(data)),
                          logical(1))

        if (!any(matches)) {
            stop(paste0("Could not determine method for ", name, " from the columns in 'data', please specify the method explicitly."))
            }

        method <- names(method_reqs)[matches][1]  # Take the first matching method
    }    

    ## Check data consistent with method:
    if(method %in% names(method_reqs)) ## This should always be true, but just in case:
        if(!(all(method_reqs[[method]] %in% names(data))))
            stop(paste0("To use '", method, "' for variable ", name,", the data must have the columns:", paste0(method_reqs[[method]],collapse = ", "), "."))   

    ##
    ## Define method functions
    ## 
    
    ## Last value carried forward
    if(method == "last"){
        data <- copy(data)[, .(id, date, value)]
        ## Avoid unlikely cases of (id,time)-duplicates
        data = data[, .(value = mean(value)), by = c("id","date")]
        out = map_grid(grid=grid,
                       data=data,
                       name=name,
                       fun_aggregate = NULL,
                       values = NULL,
                       rollforward=Inf,
                       id = id)
        return(out)
    }
    
    ## Exposed at some time in interval
    if(method == "exposed"){
        out = exposure_time(data = data,grid = gg,name = name,id = id,threshold = 0)
        return(out)
    }

    ## Exposure time
    if(method == "exposure_time"){
        out = exposure_time(data = data,grid = gg,name = name,id = id,threshold = NA)
        return(out)
    }

    ## Some event in interval least one time point in interval
    ## Should check this
    if(method == "event"){
        gg = copy(grid)
        setnames(gg,c("date","previous_date"),c("end_interval","start_interval"))
        out = discretize_timevarying_exposure(data = data,
                                              grid = gg,
                                              name = name,
                                              point_exposure = TRUE,
                                              id = id)
        return(out)
    }

    stop(paste0("Method ", method, "chosen for variable ", name, " is not implemented."))
}

