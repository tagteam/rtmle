### prepare_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (10:07) 
## Version: 
## Last-Updated: Sep 30 2024 (09:18) 
##           By: Thomas Alexander Gerds
##     Update #: 109
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
##' }
##' @return The object augmented with a new element called \code{prepared_data}.
##' @seealso rtmle_init
##' @examples
##' set.seed(112)
#' ld <- simulate_long_data(n = 101,number_epochs = 20,
#'                          beta = list(A_on_Y = -.2,
#'                          A0_on_Y = -0.3,A0_on_A = 6),
#'                                register_format = TRUE)
#' x <- rtmle_init(intervals = 8, name_id = "id",
#'                 name_outcome = "Y", name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
#' baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
#' x <- long_to_wide(x,intervals=seq(0,2000,30.45*6))
#' prepare_data(x) <- list()
#'                        
#' x$prepared_data$data
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
"prepare_data<-" <- function(x,...,value){
    ## if (length(x$protocols) == 0) stop("Object contains no protocols yet. Need at least one protocol to prepare treatment variables.")
    ## if (length(x$prepared_data) == 0 || (length(value$refit)>0 && value$refit)){
    ## x$prepared_data <- vector(length(x$protocols),mode = "list")
    ## names(x$prepared_data) <- names(x$protocols)
    ## }
    K = length(x$times)
    time_horizon = max(x$times)
    stopifnot(is.list(value))
    stopifnot(inherits(x$data$outcome_data,"data.table"))
    stopifnot(match(x$names$id,names(x$data$outcome_data),nomatch = 0)>0)
    subset = value$subset
    work_data <- x$data$outcome_data
    data.table::setkeyv(work_data,x$names$id)
    # adding the baseline covariates
    name_baseline_covariates = setdiff(names(x$data$baseline_data),x$names$id)
    if (!is.null(x$data$baseline_data)){
        stopifnot(inherits(x$data$baseline_data,"data.table"))
        stopifnot(match(x$names$id,names(x$data$baseline_data),nomatch = 0)>0)
        work_data=x$data$baseline_data[work_data,on = x$names$id]
    }
    # adding the timevarying covariates including treatment variables
    name_time_covariates <- NULL
    if (is.data.frame(x$data$timevar_data)){
        name_time_covariates <- unique(sapply(strsplit(names(x$data$timevar_data),"_"),"[",1))[-1]
        work_data <- work_data[x$data$timevar_data,on = x$names$id]
    }else{
        name_time_covariates <- names(x$data$timevar_data)
        if (length(name_time_covariates)>0){
            for (Vname in name_time_covariates){
                work_data <- work_data[x$data$timevar_data[[Vname]],on = x$names$id]
            }
        }
    }
    #
    # exclude subjects with outcome/death/censored at time zero
    # 
    Y_0 = match(paste0(x$names$outcome,"_",0),names(work_data))
    D_0 = match(paste0(x$names$competing,"_",0),names(work_data))
    C_0 = match(paste0(x$names$censoring,"_",0),names(work_data))
    excluded <- FALSE
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
    # sorting the variables for LTMLE
    work_data = work_data[,c(x$names$id, intersect(c(name_baseline_covariates,unlist(sapply(x$times, function(timepoint){
        if(timepoint == 0){
            paste0(name_time_covariates,"_",timepoint)
        } else{
            if(timepoint != x$times[K]){
                paste0(c(x$names$censoring, x$names$outcome, x$names$competing, name_time_covariates),"_",timepoint)
            } else {
                paste0(c(x$names$censoring, x$names$outcome),"_",timepoint)
            }
        }
    }))), names(work_data))), with = FALSE]
    # subset data
    if(length(subset)>0){
        subset_dt = data.table(ID = subset)
        setnames(subset_dt,"ID",x$names$id)
        work_data = work_data$data[subset_dt,on = x$names$id]
    }
    # label the variables that are constant in the subset data
    same = sapply(work_data, function(x){length(unique(x))==1})
    if(sum(same)>0){
        constant_variables <- names(work_data)[same]
    } else{
        constant_variables <- NULL}
    name_baseline_covariates <- intersect(name_baseline_covariates,names(work_data))
    ## check
    if(length(x$names$censoring)>0){
        for(col in sapply(x$times[-1], function(timepoint){paste0(x$names$censoring,"_",timepoint)})){
            set(work_data, j = col, value=as.factor(ifelse(work_data[[col]]%in%x$censored_label,"censored","uncensored")))
        }
    }
    ## Manipulation of the event nodes
    Y_nodes = unlist(lapply(x$times[-1], function(time){paste0(x$names$outcome, "_", time)}))
    D_nodes = unlist(lapply(x$times[-c(1,K)], function(time){paste0(x$names$competing, "_", time)}))
    C_nodes = unlist(lapply(x$times[-1], function(time){paste0(x$names$censoring, "_", time)}))
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
        ## The last Y_node is the last node and hence no action needed
        for(k in Y_nodes_position[-(K-1)]){
            later_nodes=setdiff((k+1):NCOL(work_data),Y_nodes_position)
            later_Y_nodes=intersect((k+1):Y_nodes_position[length(Y_nodes_position)],Y_nodes_position)
            if(any(has_outcome <- (work_data[[k]]%in%1))){
                for(l in later_nodes) {set(work_data,j=l,i=which(has_outcome),value=NA)}
                for(l in later_Y_nodes) {set(work_data,j=l,i=which(has_outcome),value=1)}
            }
        }
        if(length(x$names$competing)>0){
            for(k in D_nodes_position){
                later_nodes=setdiff((k+1):NCOL(work_data),Y_nodes_position)
                # Later outcome event nodes are set to 0
                later_Y_nodes=intersect((k+1):NCOL(work_data),Y_nodes_position)
                if(any(has_died <- (work_data[[k]]%in%1))){
                    for(l in later_nodes) {set(work_data,j=l,i=which(has_died),value=NA)}
                    for(l in later_Y_nodes) {set(work_data,j=l,i=which(has_died),value=0)}
                }
            }
        }
        ## All nodes should be NA as soon as censoring has occurred
        if(length(x$names$censoring)>0){
            for(k in C_nodes_position){
                later_nodes=(k+1):NCOL(work_data)
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
        if(length(x$names$censoring)>0){
            for(k in C_nodes_position){
                later_nodes=(k+1):NCOL(work_data)
                if(any(has_censored <- (work_data[[k]]%in%"censored"))){
                    for(l in later_nodes) {set(work_data,j=l,i=which(has_censored),value=NA)}
                }
            }
        }
    }
    L_nodes <- c(sapply(x$times, function(k) {paste0(c(name_time_covariates), "_", k)}))
    L_nodes <- L_nodes[match(L_nodes, names(work_data),nomatch = 0)!=0]
    if(length(x$names$censoring)==0) {C_nodes = NULL}
    x$prepared_data$data <- work_data[]
    x$prepared_data$name_time_covariates <- name_time_covariates
    x$prepared_data$name_baseline_covariates <- name_baseline_covariates
    x$prepared_data$name_constant_variables <- constant_variables
    x$prepared_data$Cnodes = intersect(C_nodes, names(work_data))
    x$prepared_data$Dnodes = intersect(D_nodes, names(work_data))
    x$prepared_data$Lnodes = intersect(L_nodes, names(work_data)) 
    x$prepared_data$Ynodes = intersect(Y_nodes, names(work_data))
    x
}
######################################################################
### prepare_data.R ends here
