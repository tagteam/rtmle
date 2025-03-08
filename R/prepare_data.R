### prepare_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (10:07) 
## Version: 
## Last-Updated: Mar  7 2025 (18:25) 
##           By: Thomas Alexander Gerds
##     Update #: 188
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
##' for the intervention variables (A_variables), the time dependent covariates (B_variables)
##' the outcome variables (outcome_variables), the competing risk variables (competing_variables) and the
##' censoring variables (censoring_variables) and adds them to the object. 
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
#' ld <- simulate_long_data(n = 11,number_visits = 20,
#'                          beta = list(A_on_Y = -.2,
#'                          A0_on_Y = -0.3,A0_on_A = 6),
#'                                register_format = TRUE)
#' x <- rtmle_init(intervals = 3, name_id = "id",
#'                 name_outcome = "Y", name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
#' add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
#' x <- long_to_wide(x,intervals=seq(0,2000,30.45*6))
#' prepare_data(x) <- list()
#' x$prepared_data
#' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
"prepare_data<-" <- function(x,...,value){
    ## if (length(x$protocols) == 0) stop("Object contains no protocols yet. Need at least one protocol to prepare treatment variables.")
    ## if (length(x$prepared_data) == 0 || (length(value$refit)>0 && value$refit)){
    ## x$prepared_data <- vector(length(x$protocols),mode = "list")
    ## names(x$prepared_data) <- names(x$protocols)
    ## }
    K = length(x$times)
    max_time_horizon = max(x$times)
    stopifnot(max_time_horizon>0)
    stopifnot(is.list(value))
    if (!inherits(x$data$outcome_data,"data.table"))
        if (inherits(x$data$outcome_data,"data.frame")|
            inherits(x$data$outcome_data,"tibble")){
            setDT(x$data$outcome_data)
        } else{
            stop("The object x$outcome_data must be a data.table, data.frame, or tibble")
        }
    if (match(x$names$id,names(x$data$outcome_data),nomatch = 0) == 0)
        stop("The object x$outcome_data does not have the variable '",x$names$id,"' as was initialized under x$names$id.\nPerhaps the subject identifying variable exists under a different name?")
    ## stopifnot(match(x$names$id,names(x$data$outcome_data),nomatch = 0)>0)
    #
    # FIXME: subset at the very end, that is, not in this function
    # or subset each element (outcome, baseline, timevar) before merging?
    #
    ## subset = value$subset
    # subset data
    ## if(length(subset)>0){
    ## subset_dt = data.table(ID = subset)
    ## setnames(subset_dt,"ID",x$names$id)
    ## work_data = work_data$data[subset_dt,on = x$names$id]
    ## }else{
    work_data <- x$data$outcome_data
    ##
    ## analysis of the outcome, competing and censoring variables
    ##
    outcome_variables = paste0(x$names$outcome, "_", 1:max_time_horizon)
    if(length(x$names$competing)>0){    
        competing_variables = paste0(x$names$competing, "_", 1:(max_time_horizon-1))
    }else{
        competing_variables <- NULL
    }
    if(length(x$names$censoring)>0){
        censoring_variables = paste0(x$names$censoring, "_", 1:max_time_horizon)
    }else{
        censoring_variables <- NULL
    }
    #
    # exclude subjects with outcome/death/censored events at time zero
    #
    for (y in c("outcome","competing","censoring")){
        if (match(Y_0 <- paste0(x$names[[y]],"_",0),names(work_data),nomatch = 0)>0){
            if (length(unique(work_data[[Y_0]])) > 1){
                stop(paste0("The object 'x$data$outcome_data' contains a variable at time zero called ",
                            paste0(x$names$outcome,"_",0),
                            " which is not constant. But, this variable is used to define the end-of-followup: min(",x$names$outcome,", ", x$names$competing,", ",x$names$censoring,")"))
            }else{
                work_data[[paste0(x$names[[y]],"_",0)]] <- NULL
            }
        }
    }
    data.table::setkeyv(work_data,x$names$id)
    # merging with the baseline covariates
    if (length(x$data$baseline_data)>0){
        name_baseline_covariates = setdiff(names(x$data$baseline_data),x$names$id)
        if (!is.null(x$data$baseline_data)){
            stopifnot(inherits(x$data$baseline_data,"data.table"))
            stopifnot(match(x$names$id,names(x$data$baseline_data),nomatch = 0)>0)
            work_data=x$data$baseline_data[work_data,on = x$names$id]
        }
    }else{
        name_baseline_covariates <- NULL
    }
    # merging with the timevarying covariates including treatment variables
    name_time_covariates <- NULL
    if (length(x$data$timevar_data)>0){
        if (is.data.frame(x$data$timevar_data)){
            # a time varying covariate is recognized when
            # there is a _1 version of it. That is a variable called hba_1 (root: hba) and also
            # a variable hba1c_level_1_1 (root: hba1c_level_1)
            # FIXME: there could be root variables which have a time _3 version but not a time _1 version
            name_time_covariates <- unlist(lapply(grep("_1$",names(x$data$timevar_data),value=TRUE),
                                                  function(x){substring(x,0,nchar(x)-2)}))
            work_data <- work_data[x$data$timevar_data,on = x$names$id]
        }else{
            name_time_covariates <- names(x$data$timevar_data)
            if (length(name_time_covariates)>0){
                for (Vname in name_time_covariates){
                    work_data <- work_data[x$data$timevar_data[[Vname]],on = x$names$id]
                }
            }
        }
    } else {
        name_time_covariates <- NULL
    }
    # sorting of variables should not be necessary because the formulas define the dependencies
    # but sorting is still be convenient for data inspections
    # and also our data manipulation below where we set values to NA after censoring still depends on the order 
    work_data <- work_data[,c(x$names$id, intersect(c(name_baseline_covariates,unlist(sapply(x$times, function(timepoint){
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
    # now the number of columns are final and we can search for the outcome variables 
    outcome_variables_position = match(outcome_variables, names(work_data))
    if (any(this_out <- is.na(outcome_variables_position))){
        stop(paste0("Cannot find outcome variable(s):\n",paste0(outcome_variables[this_out],collapse = ", "),"\n in x$data$outcome_data"))
    }
    if(length(x$names$censoring)>0){
        censoring_variables_position = match(censoring_variables, names(work_data))
        if (any(this_out <- is.na(censoring_variables_position))){
            stop(paste0("Cannot find censoring variable(s):\n",paste0(censoring_variables[this_out],collapse = ", "),"\n in x$data$outcome_data"))
        }
    }
    if(length(x$names$competing)>0){    
        competing_variables_position = match(competing_variables, names(work_data))
        if (any(this_out <- is.na(competing_variables_position))){
            stop(paste0("Cannot find competing risk variable(s):\n",paste0(competing_variables[this_out],collapse = ", "),"\n in x$data$outcome_data"))
        }
    }
    # label the variables that are constant in the (subset) data
    same <- sapply(work_data, function(x){length(unique(x))==1})
    if(sum(same)>0){
        constant_variables <- names(work_data)[same]
    } else{
        constant_variables <- NULL}
    name_baseline_covariates <- intersect(name_baseline_covariates,names(work_data))
    # FIXME: remove or elaborate this sanity check 
    stopifnot(nrow(work_data)>0)
    ## make sure that the order of the censoring variables are factors with ordered labels 
    ## such that we estimate the probability of being uncensored and not the probability being censored
    if(length(x$names$censoring)>0){
        for(col in censoring_variables){
            set(work_data,
                j = col,
                value=factor(work_data[[col]], levels = c(x$names$censored_label,x$names$uncensored_label)))
        }
    }
    ##
    ## due to the discretization of time it may happen that a person who has an event or dies in the interval
    ## appears also as censored in the interval. in this case we change the label of the censoring variable
    ## to uncensored 
    ##
    if (length(x$names$censoring)>0){
        for(j in 1:(max_time_horizon)){
            if((j<max_time_horizon) & !(length(x$names$competing_variables) == 0)){
                has_outcome_and_censored <- (((work_data[[outcome_variables[[j]]]]%in%1)
                    |(work_data[[competing_variables[[j]]]]%in%1))
                    & (work_data[[censoring_variables[[j]]]]%in%x$names$censored_label))
            } else{
                has_outcome_and_censored <- ((work_data[[outcome_variables_position[[j]]]]%in%1)
                    &(work_data[[censoring_variables_position[[j]]]]%in%x$names$censored_label))
            }
            if(any(has_outcome_and_censored)){
                set(work_data,j=censoring_variables[[j]],
                    i=which(has_outcome_and_censored),
                    value=x$names$uncensored_label)
            }
        }
    }
    ## All variables (except outcome and competing risk) should be NA after an event (outcome or death)
    ## outcome should be 1 and competing risk should be 0 after an outcome event
    ## outcome should be 0 and competing risk should be 1 after a competing event 
    if(max_time_horizon!= 1){
        ## The last outcome_variable is the very last variable and hence no action needed
        for(k in outcome_variables_position[-(K-1)]){
            later_variables=setdiff((k+1):NCOL(work_data),outcome_variables_position)
            later_outcome_variables=intersect((k+1):outcome_variables_position[length(outcome_variables_position)],outcome_variables_position)
            if(length(x$names$competing)>0){
                later_competing_variables=intersect((k+1):NCOL(work_data),competing_variables_position)
            }
            if(any(has_outcome <- (work_data[[k]]%in%1))){
                for(l in later_variables) {set(work_data,j=l,i=which(has_outcome),value=NA)}
                for(l in later_outcome_variables) {set(work_data,j=l,i=which(has_outcome),value=1)}
                if(length(x$names$competing)>0){
                    for(l in later_competing_variables) {set(work_data,j=l,i=which(has_outcome),value=0)}
                }
            }
        }
        if(length(x$names$competing)>0){
            for(k in competing_variables_position){
                later_variables=setdiff((k+1):NCOL(work_data),outcome_variables_position)
                later_outcome_variables=intersect((k+1):NCOL(work_data),outcome_variables_position)
                later_competing_variables=intersect((k+1):NCOL(work_data),competing_variables_position)
                if(any(has_died <- (work_data[[k]]%in%1))){
                    for(l in later_variables) {set(work_data,j=l,i=which(has_died),value=NA)}
                    for(l in later_outcome_variables) {set(work_data,j=l,i=which(has_died),value=0)}
                    for(l in later_competing_variables) {set(work_data,j=l,i=which(has_died),value=1)}
                }
            }
        }
        ## All variables should be NA as soon as censoring has occurred
        if(length(x$names$censoring)>0){
            for(k in censoring_variables_position){
                later_variables=(k+1):NCOL(work_data)
                if(any(has_censored <- (work_data[[k]]%in%x$names$censored_label))){
                    for(l in later_variables) {set(work_data,j=l,i=which(has_censored),value=NA)}
                }
            }
        }
    }else{
        # 
        # max_time_horizon = 1 we set the outcome to NA in case of censored
        #                  and same for competing risks
        #
        if(length(x$names$censoring)>0){
            for(k in censoring_variables_position){
                later_variables=(k+1):NCOL(work_data)
                if(any(has_censored <- (work_data[[k]]%in%x$names$censored_label))){
                    for(l in later_variables) {set(work_data,j=l,i=which(has_censored),value=NA)}
                }
            }
        }
    }
    ## 
    ## persons are still followed ("last_interval") if they are uncensored and free of outcome and free of competing events
    ## recall that K is the number of time points where the first time point is always 0
    ##  
    x$followup <- work_data[,c(x$names$id),with = FALSE]
    set(x$followup,j = "last_interval",value = numeric(NROW(x$followup)))
    ## if (length(censoring_variables)>0){
    ## x$followup[,uncensored := numeric(.N)]
    ## }
    if (max_time_horizon>1){
        for (j in 0:(max_time_horizon-2)){
            vital <- 1*(work_data[[outcome_variables[[j+1]]]] %in% "0")
            if (length(censoring_variables)>0){
                vital <- vital * 1*(work_data[[censoring_variables[[j+1]]]] %in% x$names$uncensored_label)
                ## UC <- 1*(work_data[[censoring_variables[[j+1]]]] %in% x$names$uncensored_label)
                ## vital <- vital * UC
                ## x$followup[,uncensored := uncensored+UC]
            }
            if (length(competing_variables)>0)
                vital <- vital * (work_data[[competing_variables[[j+1]]]] %in% "0")
            set(x$followup,j = "last_interval",value = x$followup[["last_interval"]]+vital)
        }
    }
    x$prepared_data <- work_data[]
    x$names$name_time_covariates <- name_time_covariates
    x$names$name_baseline_covariates <- name_baseline_covariates
    x$names$name_constant_variables <- constant_variables
    x
}
######################################################################
### prepare_data.R ends here
