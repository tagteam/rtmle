### prepare_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (10:07) 
## Version: 
## Last-Updated: Jul 29 2024 (09:57) 
##           By: Thomas Alexander Gerds
##     Update #: 47
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' If the object contains data in long format these are
##' transformed to wide format. The wide format data are sorted and
##' checked for the analysis and added as \code{x$prepared_data}.
##'
##' As a wanted side effect the function identifies the variables
##' for the intervention nodes (A_nodes), the time dependent covariates (B_nodes)
##' the outcome nodes (Y_nodes), the competing risk nodes (D_nodes) and the
##' censoring nodes (C_nodes) and adds them to the object. 
##' @title Preparing a targeted minimum loss based analysis
##' @param x Object which contains long format or wide format data
##' @param ... not used
##' @param value A list with the following forced elements:
##' \itemize{
##' \item \code{"intervals"} The time intervals
##' \item \code{"treatment_variables"} The treatment variables
##' }
##' @return The object augmented with a new element called \code{prepared_data}.
##' @seealso rtmle_init
##' @examples
##' set.seed(112)
#' ld <- simulate_long_data(n = 100,number_epochs = 20,
#'                          beta = list(A_on_Y = -.2,
#'                          A0_on_Y = -0.3,A0_on_A = 6),
#'                                register_format = TRUE)
#' x <- rtmle_init(intervals = 8, name_id = "id",
#'                 name_outcome = "Y", name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' add_long_data(x) <- ld
#' prepare_data(x) <- list(treatment_variables = "A",
#'                        reset = TRUE,
#'                        intervals = seq(0,2000,30.45*6))
#' x$prepared_data
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
"prepare_data<-" <- function(x,...,value){
    stopifnot(is.list(value))
    stopifnot(all(c("intervals","treatment_variables")%in%names(value)))
    intervals = value$intervals
    reset = value$reset
    if (length(reset) == 0) reset = FALSE
    subset = value$subset
    treatment_variables = value$treatment_variables
    time_horizon = max(x$times)
    K = length(x$times)
    if (length(x$data) == 0 || reset){
        if (length(x$long_data) == 0) {stop("Object contains no data.")}
        # FIXME: this next line does not copy the data, right?
        x$data$baseline_data = x$long_data$baseline_data
        if (missing(intervals)) {stop("Need time intervals to map the long data.")}
        pop <- x$long_data$baseline_data[,x$name_id,with = FALSE]
        pop[,start_followup := rep(0,.N)]
        if (!("competingrisk_date" %in% names(x$long_data$competingrisk_data)))
            setnames(x$long_data$competingrisk_data,"date","competingrisk_date")
        pop <- x$long_data$competingrisk_data[pop,on = x$name_id]
        if (!("censored_date" %in% names(x$long_data$censored_data)))
            setnames(x$long_data$censored_data,"date","censored_date")
        pop <- x$long_data$censored_data[pop,on = x$name_id]
        if (!("outcome_date" %in% names(x$long_data$outcome_data)))
            setnames(x$long_data$outcome_data,"date","outcome_date")
        pop <- x$long_data$outcome_data[pop,on = x$name_id]
        pop[,end_followup := pmin(censored_date,competingrisk_date,outcome_date,na.rm = TRUE)]
        if (any(is.na(pop$end_followup)))stop("Missing values in end of followup information")
        grid <- pop[,.(date=start_followup+intervals, end = end_followup),by=eval(as.character(x$name_id))]
        grid[,interval:=0:(length(intervals)-1),by=eval(as.character(x$name_id))]
        grid <- pop[,.SD,.SDcols = c(x$name_id)][grid,on = eval(as.character(x$name_id))]
        length_interval=unique(round(diff(intervals),0))
        # now awkwardly reset the names
        setnames(x$long_data$competingrisk_data,"competingrisk_date","date")
        setnames(x$long_data$censored_data,"censored_date","date")
        setnames(x$long_data$outcome_data,"outcome_date","date")
        # FIXME: allow for data after the end of followup?
        ## grid <- grid[date<=end+length_interval]
        x$data$outcome_data <- widen_outcome(outcome_name = x$name_outcome,
                                             outcome_data = x$long_data$outcome_data,
                                             competingrisk_data = x$long_data$competingrisk_data,
                                             censored_data = x$long_data$censored_data,
                                             grid = grid,
                                             fun.aggregate = NULL,
                                             id = x$name_id)
        # FIXME: if more than one instance (such as treatment) happens
        #        within one interval
        #        an aggregate function determines the value inside
        #        the interval
        funA <- function(x){1*(sum(x)>0)}
        x$data$treatment_data <- map_grid(grid=grid,data=x$long_data$treatment_data,name=treatment_variables,fun.aggregate = funA,rollforward=(length_interval - 1),id = x$name_id)
        # FIXME: names of timevar data are not just L
        x$data$timevar_data <- map_grid(grid=grid,data=x$long_data$timevar_data,name="L",fun.aggregate = funA,rollforward=(length_interval - 1),id = x$name_id)
    } else{
        stopifnot(inherits(x$data$treatment_data,"data.table"))
        stopifnot(inherits(x$data$outcome_data,"data.table"))
        stopifnot(match(x$name_id,names(x$data$treatment_data),nomatch = 0)>0)
        stopifnot(match(x$name_id,names(x$data$outcome_data),nomatch = 0)>0)
    }
    data.table::setkeyv(x$data$treatment_data,x$name_id)
    data.table::setkeyv(x$data$outcome_data,x$name_id)
    work_data=x$data$outcome_data[x$data$treatment_data,on = x$name_id]
    # deal with outcome/death/censored at index
    Y_0 = match(paste0(x$name_outcome,"_",0),names(work_data))
    D_0 = match(paste0(x$name_competing,"_",0),names(work_data))
    C_0 = match(paste0(x$name_censoring,"_",0),names(work_data))
    excluded = FALSE
    if (!is.na(Y_0)){
        if(!is.na(D_0)&!is.na(C_0)){
            excluded <- (work_data[[Y_0]]%in%1)|(work_data[[D_0]]%in%1)|(work_data[[C_0]]%in%x$censored_label)
        }
        if(!is.na(D_0)){
            excluded <- (work_data[[Y_0]]%in%1)|(work_data[[D_0]]%in%1)
        }
        if(!is.na(C_0)){
            excluded <- (work_data[[Y_0]]%in%1)|(work_data[[C_0]]%in%x$censored_label)
        }
    }else{
        if(!is.na(D_0)&!is.na(C_0)){
            excluded <- (work_data[[D_0]]%in%1)|(work_data[[C_0]]%in%x$censored_label)
        }
        if(!is.na(D_0)){
            excluded <- (work_data[[D_0]]%in%1)
        }
        if(!is.na(C_0)){excluded <- (work_data[[C_0]]%in%x$censored_label)
        }
    }
    if (any(excluded)){
        work_data = work_data[!excluded]
        stopifnot(nrow(work_data)>0)
    }
    # adding the baseline covariates
    if (!is.null(x$data$baseline_data)){
        stopifnot(inherits(x$data$baseline_data,"data.table"))
        stopifnot(match(x$name_id,names(x$data$baseline_data),nomatch = 0)>0)
        work_data=x$data$baseline_data[work_data,on = x$name_id]
    }
    # add time covariates
    # first remove outcome if overlap
    if (length(x$data$timevar_data)>0){
        stopifnot(inherits(x$data$timevar_data,"data.table"))
        stopifnot(match(x$name_id,names(x$data$timevar_data),nomatch = 0)>0)
        if (length((outcome_overlap <- grep(paste0(x$name_outcome,"_"),names(x$data$timevar_data)))>0)){
            timevar_data <- x$data$timevar_data[,-outcome_overlap, with=FALSE]}
        else
            timevar_data = x$data$timevar_data
        data.table::setkeyv(timevar_data,x$name_id)
        work_data=timevar_data[work_data, on = x$name_id]
        name_time_covariates = unlist(lapply(grep("_0",names(timevar_data),value=TRUE),
                                             function(x){substring(x,0,nchar(x)-2)}))
    }else{
        name_time_covariates <- NULL
    }
    name_baseline_covariates = setdiff(names(x$data$baseline_data),x$name_id)
    # sorting the variables for LTMLE
    work_data = work_data[,c(x$name_id, intersect(c(name_baseline_covariates,unlist(sapply(x$times, function(timepoint){
        if(timepoint == 0){
            paste0(c(name_time_covariates, treatment_variables),"_",timepoint)
        } else{
            if(timepoint != x$times[K]){
                paste0(c(x$name_censoring, x$name_outcome, x$name_competing, name_time_covariates, treatment_variables),"_",timepoint)
            } else {
                paste0(c(x$name_censoring, x$name_outcome),"_",timepoint)
            }
        }
    }))), names(work_data))), with = FALSE]
    # subset data
    if(length(subset)>0){
        subset_dt = data.table(ID = subset)
        setnames(subset_dt,"ID",x$name_id)
        work_data = work_data$data[subset_dt,on = x$name_id]
    }
    # label the variables that are constant in the subset data
    same = sapply(work_data, function(x){length(unique(x))==1})
    if(sum(same)>0){
        constant_variables <- names(work_data)[same]
    } else{
        constant_variables <- NULL}
    name_baseline_covariates <- intersect(name_baseline_covariates,names(work_data))
    ## check
    if(length(x$name_censoring)>0){
        for(col in sapply(x$times[-1], function(timepoint){paste0(x$name_censoring,"_",timepoint)})){
            set(work_data, j = col, value=as.factor(ifelse(work_data[[col]]%in%x$censored_label,"censored","uncensored")))
        }
    }
    ## Manipulation of the event nodes
    A_nodes = unlist(lapply(x$times[-K], function(time){paste0(treatment_variables, "_", time)}))
    Y_nodes = unlist(lapply(x$times[-1], function(time){paste0(x$name_outcome, "_", time)}))
    D_nodes = unlist(lapply(x$times[-c(1,K)], function(time){paste0(x$name_competing, "_", time)}))
    C_nodes = unlist(lapply(x$times[-1], function(time){paste0(x$name_censoring, "_", time)}))
    # if A,B then B_0 is obsolete because A0=1-B0
    if (length(treatment_variables) == 2)
        A_nodes <- A_nodes[A_nodes != paste0(treatment_variables[[2]],"_0")]
    A_nodes_position = match(A_nodes, names(work_data))
    Y_nodes_position = match(Y_nodes, names(work_data))
    D_nodes_position = match(D_nodes, names(work_data))
    C_nodes_position = match(C_nodes, names(work_data))
    ## Adjust data depending on censoring/event/competing event with NA
    for(q in 1:(K-1)){
        if(q<(K-1)){
            has_outcome_or_death_and_censored = (((work_data[[Y_nodes_position[[q]]]]%in%"1")|(work_data[[D_nodes_position[[q]]]]%in%"1"))&
                                                 (work_data[[C_nodes_position[[q]]]]%in%"censored"))
        } else{
            has_outcome_or_death_and_censored = ((work_data[[Y_nodes_position[[q]]]]%in%1)&(work_data[[C_nodes_position[[q]]]]%in%"censored"))
        }
        if(any(has_outcome_or_death_and_censored)){
            set(work_data,j=C_nodes_position[[q]],i=which(has_outcome_or_death_and_censored),value="uncensored")
        }
    }
    ## All nodes (except outcome and competing risk) should be NA after an event (outcome or death)
    if(time_horizon!= 1){
        for(k in Y_nodes_position[-(K-1)]){
            later_nodes=setdiff((k+1):NCOL(work_data),Y_nodes_position)
            later_Y_nodes=intersect((k+1):Y_nodes_position[length(Y_nodes_position)],Y_nodes_position)
            # later_C_nodes=intersect((k+1):NCOL(work_data),C_nodes_position)
            # later_D_nodes=intersect((k+1):NCOL(work_data),D_nodes_position)
            if(any(has_outcome <- (work_data[[k]]%in%1))){
                for(l in later_nodes) {set(work_data,j=l,i=which(has_outcome),value=NA)}
                for(l in later_Y_nodes) {set(work_data,j=l,i=which(has_outcome),value=1)}
                # for(l in later_C_nodes) {set(work_data,j=l,i=which(has_outcome),value="uncensored")}
                # for(l in later_D_nodes) {set(work_data,j=l,i=which(has_outcome),value=0)}
            }
        }
        if(length(x$name_competing)>0){
            for(k in D_nodes_position){
                later_nodes=setdiff((k+1):NCOL(work_data),Y_nodes_position)
                # Later outcome event nodes are set to 0
                later_Y_nodes=intersect((k+1):NCOL(work_data),Y_nodes_position)
                # later_C_nodes=intersect((k+1):NCOL(work_data),C_nodes_position)
                # later_D_nodes=intersect((k+1):NCOL(work_data),D_nodes_position)
                if(any(has_died <- (work_data[[k]]%in%1))){
                    for(l in later_nodes) {set(work_data,j=l,i=which(has_died),value=NA)}
                    for(l in later_Y_nodes) {set(work_data,j=l,i=which(has_died),value=0)}
                    # for(l in later_C_nodes) {set(work_data,j=l,i=which(has_died),value="uncensored")}
                    # for(l in later_D_nodes) {set(work_data,j=l,i=which(has_died),value=1)}
                }
            }
        }
        ## All nodes should be NA as soon as censoring has occurred
        if(length(x$name_censoring)>0){
            for(k in C_nodes_position){
                later_nodes=(k+1):NCOL(work_data)
                ## print(k)
                ## print(names(work_data)[later_nodes])
                if(any(has_censored <- (work_data[[k]]%in%"censored"))){
                    for(l in later_nodes) {set(work_data,j=l,i=which(has_censored),value=NA)}
                }
            }
        }
    }else{
        # 
        # time_horizon = 1 we set the outcome to NA in case of censored
        #                  and same for competing risks
        #
        if(length(x$name_censoring)>0){
            for(k in C_nodes_position){
                later_nodes=(k+1):NCOL(work_data)
                ## print(names(work_data)[later_nodes])
                if(any(has_censored <- (work_data[[k]]%in%"censored"))){
                    for(l in later_nodes) {set(work_data,j=l,i=which(has_censored),value=NA)}
                }
            }
        }
    }
    L_nodes <- c(sapply(x$times, function(k) {paste0(c(name_time_covariates), "_", k)}))
    L_nodes <- L_nodes[match(L_nodes, names(work_data),nomatch = 0)!=0]
    if(length(x$name_censoring)==0) {C_nodes = NULL}
    x$prepared_data <- work_data[]
    x$name_time_covariates <- name_time_covariates
    x$name_baseline_covariates <- name_baseline_covariates
    x$name_constant_variables <- constant_variables
    x$Anodes = intersect(A_nodes, names(work_data))
    x$Cnodes = intersect(C_nodes, names(work_data))
    x$Dnodes = intersect(D_nodes, names(work_data))
    x$Lnodes = intersect(L_nodes, names(work_data)) 
    x$Ynodes = intersect(Y_nodes, names(work_data))
    x
}
######################################################################
### prepare_data.R ends here
