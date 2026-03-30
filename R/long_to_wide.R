### long_to_wide.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 22 2024 (14:07) 
## Version: 
## Last-Updated: mar 30 2026 (14:08) 
##           By: Thomas Alexander Gerds
##     Update #: 287
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Transform data from long to wide format
#'
#' This function is used to prepare a discretized time analysis of longituginal data.
#' 
#' The start_followup_date is the calendar date where the follow up
#' starts or can be zero if the followup data are readily on the time
#' on study scale
#' @param x object of class \code{rtmle}
#' @param start_followup_date Date where the followup starts. This date is substracted from the dates of
#' the outcomes, competing risks, censoring (end-of-followup), and dates where the time-varying covariates change.
#' If missing it is assumed to be zero for all subjects in which case the dates of the outcomes,
#' competing risks, censoring (end-of-followup), and time-varying covariates must be
#' given on the time on study scale.
#' @param mappings Named list of instructions for how to map the continuous time history of a variable in long-format
#' onto the discrete time scale. The names must match \code{names(x$long_data$timevar_data)}. 
#' @param verbose Logical. If \code{FALSE} suppress all messages.
#' @param ... alternative way to specify elements of \code{mappings}.
#' @details The function discretizes dates of events and concomittant marker information.
#' Multiple wide format variables may result from a single long-format variable.
#' @examples
#' set.seed(17)
#' x <- rtmle_init(time_grid = 0:2,name_id = "id",
#'                 name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- add_long_data(x,outcome_data=data.frame(id=c(1,3,5),date=0.1,1.3,1.9),
#'                      censored_data=data.frame(id=c(2,6),date=c(0.5,1.5)),
#'                     competing_data=data.frame(id=c(4,7),date=c(1.1,1.8)),
#'                     timevar_data=list(A=data.frame(id=c(1,1,1,2,2,3,4,4),
#'                                            start_date=c(0,.3,1.3,0,1,0,1.1,1.5),
#'                                              end_date=c(.3,.6,2,1,2,.5,1.4,1.9)),
#'                                       V=data.frame(id=c(1,1,1,2,2,3,4,4),
#'                                            start_date=c(0,.5,1,0,1.3,0,0,0.4),
#'                                              end_date=c(.3,.8,1.8,.5,1.4,2,.4,1.6)),
#'                                       L=data.frame(id=c(1,2,2,2,3,4,5,6,7,7),
#'                                                  date=c(0,0,.25,.75,0,0,0,0,0,1.4),
#'                                                 value=c(4,35,27.7,28.2,8.8,2,3.1,7,7.7,8.4))))
#' x <- add_baseline_data(x,data=data.frame(id=1:7,age=40:46))
#' x <- long_to_wide(x,L=list(method="locf"))
#' x$data$timevar_data$L
#' x$data$timevar_data$V
#' x <- long_to_wide(x,L=list(method="measurement",fun_aggregate="median"))
#' x$data
#' 
#' @export
long_to_wide <- function(x,
                         start_followup_date,
                         mappings,
                         verbose = TRUE,
                         ...){
    previous_date = interval = end_followup = censored_date =  competing_date = outcome_date = NULL
    breaks <- x$time_grid_labels
    # check if baseline_data exists
    if (length(x$data$baseline_data) == 0){
        stop("To discretize the long format data we need the baseline data stored as 'x$data$baseline_data' with the subject id variable.\nUse the function 'add_baseline_data' to add this.")
    }
    # check if long_data exists 
    if (length(x$long_data) == 0) {
        stop("x$long_data has length 0.")
    }
    Vnames <- names(x$long_data$timevar_data)
    if (length(Vnames)>0){
        # FIXME: this should have been taken care of by the add_long_data function 
        if (any(duplicated(Vnames))) stop("Duplicated names found in names(x$long_data$timevar_data). Variables must have distinct names.")
        #
        # resolve long-to-wide mappings
        # 
        if(missing(mappings)) mappings <- NULL
        ## Add fun supplied by dots
        dots <- list(...)
        if ((hit <- match("breaks",names(dots),nomatch = 0))>0){
            warning("rtmle::long_to_wide: Argument 'breaks' is obsolete." )
            dots <- dots[-hit]
        }
        if (length(dots) == 0) {
            dots <- NULL
        } else {
            bad <- is.null(names(dots)) || any(names(dots) == "")
            if (bad) stop("All ... arguments must be named.")
        }
        mappings <- c(mappings,dots)
        if (length(unused <- setdiff(names(mappings),Vnames))>0){
            warning("The following functions that should map from long to wide are not applicable because they do not occur in x$long_data:\n",
                    paste0(unused,collapse = ", "))   
        }
        ## Known methods for mapping long to wide data
        long_to_wide_mappings <- list(list(method = "measurement",fun = "discretize",columns = c("date","value"),args = list(lookback_window = Inf)),
                                      list(method = "locf",fun = "discretize",columns = c("date","value"),args = list(lookback_window = Inf)),
                                      list(method = "event",fun = "discretize",columns = "date"),
                                      list(method = "event_interval",fun = "discretize",columns = "date"),
                                      list(method = "any_exposure",fun = "discretize",columns = c("start_date","end_date"), args = list(threshold = 0)),
                                      list(method = "exposure_time",fun = "discretize",columns = c("start_date","end_date"), args = list(threshold = 0.5,relative = FALSE)),
                                      list(method = "exposure_percent",fun = "discretize",columns = c("start_date","end_date"), args = list(threshold = 0.5,relative = TRUE)))
        known_methods <- sapply(long_to_wide_mappings,function(x)x$method)
        names(long_to_wide_mappings) <- known_methods
        for (Vname in intersect(names(mappings),Vnames)){
            if (is.character(mappings[[Vname]]) && (mappings[[Vname]] %chin% known_methods)){
                mappings[[Vname]] <- long_to_wide_mappings[[mappings[[Vname]]]]
            }else{
                V_method <- mappings[[Vname]]$method
                if (is.null(V_method)) stop("You need to specify at least one method for how to map variable ",Vname," to wide format")
                ## combine arguments with default arguments
                if (V_method %chin% known_methods){
                    mapping_Vname <- c(mappings[[Vname]],long_to_wide_mappings[[V_method]])
                    mapping_Vname <- mapping_Vname[!duplicated(names(mapping_Vname))]
                    mappings[[Vname]] <- mapping_Vname
                }else{
                    ## user-defined mapping
                    stopifnot(is.function(mappings[[Vname]]$fun))
                }
            }
            ## Check if data is consistent with mapping
            if(!(all(mappings[[Vname]]$columns %chin% names(x$long_data$timevar_data[[Vname]])))){
                stop(paste0("To apply the long-to-wide mapping method '",mappings[[Vname]]$method,"' on variable ",Vname,",\nthe data for ",Vname," provided as an element of x$long_data must have the following columns:\n",paste0(mappings[[Vname]]$columns,collapse = ", "),"."))
            }
        }
        ## for variables that do not occur in mappings or ... we try to find a suitable mapping
        if (length(unmapped <- setdiff(Vnames,names(mappings)))>0){
            unmapped_methods <- lapply(unmapped, function(v) {
                ## Attempt to guess mapping if not specified:
                guessed_mapping <- NULL
                for (m in long_to_wide_mappings) {
                    if (all(m$columns %in% names(x$long_data$timevar_data[[v]]))) {
                        guessed_mapping <- m
                        break
                    }
                }
                if (is.null(guessed_mapping)) {
                    stop(paste0("Could not determine method for timevarying variable ", v, " from the columns in 'data', please specify the method explicitly."))
                }
                guessed_mapping
            })
            names(unmapped_methods) <- unmapped
            mappings <- c(mappings,unmapped_methods)
        }
        #
        # check if all dates have the same format
        #
        tv_date_formats <- sapply(Vnames,function(v){
            if ("date"%in%names(x$long_data$timevar_data[[v]])){
                class(x$long_data$timevar_data[[v]][["date"]])
            }else{
                uclass <- unique(c(class(x$long_data$timevar_data[[v]][["start_date"]]),
                                   class(x$long_data$timevar_data[[v]][["end_date"]])))
                stopifnot(length(uclass) == 1)
                uclass
            }
        })
    }else{
        tv_date_formats <- NULL
    }
    ocd_date_formats <- sapply(intersect(c("outcome_data","censored_data","competing_data"),names(x$long_data)),function(v){
        class(x$long_data[[v]][["date"]])
    })
    if (length(unique(c(tv_date_formats,ocd_date_formats))) > 1){
        ctab <- data.table(
            "variable" = names(c(tv_date_formats,ocd_date_formats)),
            "storage" = rep(c("long_data$timevar_data","long_data"),c(length(tv_date_formats),length(ocd_date_formats))),
            "class_of_date_variable" = c(tv_date_formats,ocd_date_formats))
        cat("\n\nAnalysis of date variables:\n\n")
        print(ctab)
        cat("\n")
        stop(paste0("All date variables in long format data need to have the same class (either numeric or Date)."))
    }
    #
    # resolve time zero
    # 
    if (missing(start_followup_date)||start_followup_date[1] == 0){
        if (missing(start_followup_date)){
            x$details <- c(x$details,list("Missing start_followup_date variable, for now assume 0, which implies that all event times must be given on the scale: 'time since 0'."))
        }
        pop <- x$data$baseline_data[,c(x$names$id),with = FALSE][,start_followup_date := rep(0,.N)]
    }else{
        if (!is.character(start_followup_date) || match(start_followup_date,names(x$data$baseline_data),nomatch = 0) == 0){
            stop("Argument start_followup_date must be the name (as character) of a variable in x$data$baseline_data")
        }
        x$names$start_followup_date <- start_followup_date
        pop <- x$data$baseline_data[,c(x$names$id,start_followup_date),with = FALSE]
        setnames(pop,start_followup_date,"start_followup_date")
    }
    #
    # outcome, censored and competing risk define end-of-followup
    #
    if (length(x$names$competing)>0 & length(x$long_data$competing_data)>0){
        if (length(x$long_data$outcome_data)>0){
            if (!("competing_date" %in% names(x$long_data$competing_data)))
                setnames(x$long_data$competing_data,"date","competing_date")
        }
        pop <- x$long_data$competing_data[pop,on = x$names$id]
    }else{
        pop[,competing_date := as.Date(Inf)]
    }
    if (length(x$names$censoring)>0 & length(x$long_data$censored_data)>0){
        if (!("censored_date" %in% names(x$long_data$censored_data)))
            setnames(x$long_data$censored_data,"date","censored_date")
        pop <- x$long_data$censored_data[pop,on = x$names$id]
    }else{
        pop[,censored_date := as.Date(Inf)]
    }
    if (!("outcome_date" %in% names(x$long_data$outcome_data)))
        setnames(x$long_data$outcome_data,"date","outcome_date")
    pop <- x$long_data$outcome_data[pop,on = x$names$id]
    pop[,end_followup := pmin(censored_date,competing_date,outcome_date,na.rm = TRUE)]
    pop[,c(x$names$id,"start_followup_date","end_followup"),with = FALSE]
    if (any(is.na(pop$end_followup)))stop("Missing values in end of followup information")
    id_colname <- x$names$id
    grid <- pop[,data.table::data.table(date=start_followup_date+breaks, end_followup = end_followup),by=id_colname]
    # Note: dates after the last grid time are not used
    #       FIXME: should we warn if there are such dates?
    grid[,previous_date := c(0,date[-.N]),by = id_colname]
    grid <- grid[previous_date<=end_followup]
    grid[,end_followup := NULL]
    grid[,interval:=0:(.N-1),by=id_colname]
    grid <- pop[,.SD,.SDcols = c(x$names$id)][grid,on = id_colname]
    # now awkwardly reset the names
    if (length(x$names$competing)>0 & length(x$long_data$competing_data)>0){
        setnames(x$long_data$competing_data,"competing_date","date")
    }
    if (length(x$names$censoring)>0 & length(x$long_data$censored_data)>0){
        setnames(x$long_data$censored_data,"censored_date","date")
    }
    setnames(x$long_data$outcome_data,"outcome_date","date")
    # map outcome information to discrete time scale
    x$data$outcome_data <- widen_outcome(x,grid = grid,fun_aggregate = NULL)
    ## 
    ## Time dependent variables including treatment
    ##
    ## Construct suffix to distinguish wide_format variables that are based on the same timevar_data
    setnames(grid,c("previous_date","date"),c("start_interval","end_interval"))
    data.table::setcolorder(grid,c(x$names$id,"interval","start_interval","end_interval"))
    if (length(Vnames)>0){
        for (Vname in names(mappings)){
            m <- mappings[[Vname]]
            args <- list(method = m$method,
                         data = x$long_data$timevar_data[[Vname]],
                         grid = grid,
                         name = Vname,
                         id = x$names$id,
                         lookback_window = switch(m$method,
                                                  "event" = Inf,
                                                  "locf" = Inf, NA),
                         values = c(1, 0),
                         fun_aggregate = m$fun_aggregate,
                         fill = NA)
            args <- c(args,
                      mappings[[Vname]]$args)
            args <- args[!(duplicated(names(args)))]
            x$data$timevar_data[[Vname]] <- do.call(m$fun,args)
        }
    }
    return(x)
}

######################################################################
### long_to_wide.R ends here
