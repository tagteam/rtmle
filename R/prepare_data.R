### prepare_data.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (10:07)
## Version:
## Last-Updated: sep 10 2025 (08:31) 
##           By: Thomas Alexander Gerds
##     Update #: 240
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
##' @param x Object initialized with \code{rtmle_init} which contains wide format data.
##' @param ... not used (not yet)
##' \itemize{
##' \item \code{"intervals"} The time intervals
##' }
##' @return The object augmented with a new element called \code{prepared_data}.
##' @seealso rtmle_init
##' @examples
#' set.seed(112)
#' ld <- simulate_long_data(n = 11,number_visits = 20,
#'                          beta = list(A_on_Y = -.2,
#'                          A0_on_Y = -0.3,A0_on_A = 6),
#'                                register_format = TRUE)
#' x <- rtmle_init(intervals = 3, name_id = "id",
#'                 name_outcome = "Y", name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- add_long_data(x,
#'                    outcome_data=ld$outcome_data,
#'                    censored_data=ld$censored_data,
#'                    competing_data=ld$competing_data,
#'                    timevar_data=ld$timevar_data)
#' x <- add_baseline_data(x,data=ld$baseline_data)
#' x <- long_to_wide(x,intervals=seq(0,2000,30.45*6))
#' x <- prepare_data(x)
#' x$prepared_data
#'
##' @export
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
prepare_data <- function(x,...){
    
    # check object x
    K = length(x$times)
    max_time_horizon = max(x$times)
    stopifnot(max_time_horizon>0)

    # names of outcome, competing and censoring variables
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
    if (NROW(x$prepared_data)>0){
        warning("Overwriting existing prepared data object")
    }
    # FIXME: are these checks obsolete, now that the add_data function perform these checks?
    check_id <- function(data,id,allow_null = TRUE, name = ""){
        if (NROW(data)>0){
            if (!inherits(data,"data.table"))
                if (inherits(data,"data.frame")){
                    data.table::setDT(data)
                } else{
                    warning("The object x$data$outcome_data must be a data.table, data.frame, or tibble")
                }
            if (match(id,names(data),nomatch = 0) == 0){
                stop(paste0("The object ",name," does not have the variable '",
                            id,
                            "' as was initialized under x$names$id.\n",
                            "Perhaps the subject identifying variable exists under a different name?"))
            }
        }else{
            if (!allow_null[[1]]){
                stop(paste0("Data ",name," has zero rows"))
            }
        }
    }
    check_id(data = x$data$outcome_data, id = x$names$id,allow_null = FALSE,name = "x$data$outcome_data")
    check_id(data = x$data$baseline_data, id = x$names$id,name = "x$data$baseline_data")
    # FIXME: do something else when data were added via add_wide_data
    ## nix = lapply(names(x$data$timevar_data),function(nn){
    ## check_id(data = x$data$timevar_data[[nn]], id = x$names$id,name = paste0("x$data$timevar_data$",nn))
    ## })
    ## rm(nix)
    # initialize dataset with outcome data
    prepared_data <- x$data$outcome_data
    # sort by id
    data.table::setDT(prepared_data)
    data.table::setkeyv(prepared_data,x$names$id)
    #
    # make sure that the user has excluded all subjects with outcome/death/censored events at time zero
    #
    for (y in c("outcome","competing","censoring")){
        if (match(Y_0 <- paste0(x$names[[y]],"_",0),names(prepared_data),nomatch = 0)>0){
            if (length(unique(prepared_data[[Y_0]])) > 1){
                stop(paste0("The object 'x$data$outcome_data' contains a variable at time zero called ",
                            Y_0,
                            " which is not constant. But, this variable is used to define the end-of-followup: min(",
                            x$names$outcome,
                            ", ",
                            x$names$competing,
                            ", ",
                            x$names$censoring,
                            ").\nThe current estimation framework cannot deal with outcome and end of followup events at baseline.\n"))
            }else{
                prepared_data[[paste0(x$names[[y]],"_",0)]] <- NULL
            }
        }
    }
    # merging with the baseline covariates
    if (length(x$data$baseline_data)>0){
        name_baseline_covariates = setdiff(names(x$data$baseline_data),c(x$names$id,x$names$start_followup_date))
        if (!is.null(x$data$baseline_data)){
            stopifnot(inherits(x$data$baseline_data,"data.table"))
            stopifnot(match(x$names$id,names(x$data$baseline_data),nomatch = 0)>0)
            prepared_data=x$data$baseline_data[prepared_data,on = x$names$id]
            # remove start of follow_up date because this is not a covariate
            if (length(x$names$start_followup_date)>0)
                prepared_data[[x$names$start_followup_date]] <- NULL
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
            name_time_covariates <- unlist(lapply(grep("_[1-9]+$",names(x$data$timevar_data),value=TRUE),
                                                  function(x){substring(x,0,nchar(x)-2)}))
            prepared_data <- prepared_data[x$data$timevar_data,on = x$names$id]
        }else{
            name_time_covariates <- names(x$data$timevar_data)
            if (length(name_time_covariates)>0){
                for (Vname in name_time_covariates){
                    prepared_data <- prepared_data[x$data$timevar_data[[Vname]],on = x$names$id]
                }
            }
        }
        non_variables <- setdiff(unlist(lapply(name_time_covariates,function(x)paste0(x,"_",0:K))),
                                 names(prepared_data))
    } else {
        name_time_covariates <- NULL
        non_variables <- NULL
    }
    # sorting of variables should not be necessary because the formulas define the dependencies
    # but sorting is still be convenient for data inspections
    # and also our data manipulation below where we set values to NA after censoring still depends on the order
    prepared_data <- prepared_data[,c(x$names$id, intersect(c(name_baseline_covariates,unlist(sapply(x$times, function(timepoint){
        if(timepoint == 0){
            paste0(name_time_covariates,"_",timepoint)
        } else{
            if(timepoint != x$times[K]){
                paste0(c(x$names$censoring, x$names$outcome, x$names$competing, name_time_covariates),"_",timepoint)
            } else{
                paste0(c(x$names$censoring, x$names$outcome),"_",timepoint)
            }
        }
    }))), names(prepared_data))), with = FALSE]
    # now the number of columns are final and we can search for the outcome variables
    outcome_variables_position = match(outcome_variables, names(prepared_data))
    if (any(this_out <- is.na(outcome_variables_position))){
        stop(paste0("Cannot find outcome variable(s):\n",paste0(outcome_variables[this_out],collapse = ", "),"\n in x$data$outcome_data"))
    }
    if(length(x$names$censoring)>0){
        censoring_variables_position = match(censoring_variables, names(prepared_data))
        if (any(this_out <- is.na(censoring_variables_position))){
            stop(paste0("Cannot find censoring variable(s):\n",paste0(censoring_variables[this_out],collapse = ", "),"\n in x$data$outcome_data"))
        }
    }
    if(length(x$names$competing)>0){
        competing_variables_position = match(competing_variables, names(prepared_data))
        if (any(this_out <- is.na(competing_variables_position))){
            stop(paste0("Cannot find competing risk variable(s):\n",paste0(competing_variables[this_out],collapse = ", "),"\n in x$data$outcome_data"))
        }
    }
    # label the variables that are constant in the (subset) data
    same <- sapply(prepared_data, function(x){length(unique(x))==1})
    if(sum(same)>0){
        constant_variables <- names(prepared_data)[same]
    } else{
        constant_variables <- NULL}
    name_baseline_covariates <- intersect(name_baseline_covariates,names(prepared_data))
    # FIXME: remove or elaborate this sanity check
    stopifnot(nrow(prepared_data)>0)
    ## make sure that the order of the censoring variables are factors with ordered labels
    ## such that we estimate the probability of being uncensored and not the probability being censored
    if(length(x$names$censoring)>0){
        for(col in censoring_variables){
            set(prepared_data,
                j = col,
                value=factor(prepared_data[[col]], levels = c(x$names$censored_label,x$names$uncensored_label)))
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
                has_outcome_and_censored <- (((prepared_data[[outcome_variables[[j]]]]%in%1)
                    |(prepared_data[[competing_variables[[j]]]]%in%1))
                    & (prepared_data[[censoring_variables[[j]]]]%in%x$names$censored_label))
            } else{
                has_outcome_and_censored <- ((prepared_data[[outcome_variables_position[[j]]]]%in%1)
                    &(prepared_data[[censoring_variables_position[[j]]]]%in%x$names$censored_label))
            }
            if(any(has_outcome_and_censored)){
                set(prepared_data,j=censoring_variables[[j]],
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
            later_variables=setdiff((k+1):NCOL(prepared_data),outcome_variables_position)
            later_outcome_variables=intersect((k+1):outcome_variables_position[length(outcome_variables_position)],outcome_variables_position)
            if(length(x$names$competing)>0){
                later_competing_variables=intersect((k+1):NCOL(prepared_data),competing_variables_position)
            }
            if(any(has_outcome <- (prepared_data[[k]]%in%1))){
                for(l in later_variables) {set(prepared_data,j=l,i=which(has_outcome),value=NA)}
                for(l in later_outcome_variables) {set(prepared_data,j=l,i=which(has_outcome),value=1)}
                if(length(x$names$competing)>0){
                    for(l in later_competing_variables) {set(prepared_data,j=l,i=which(has_outcome),value=0)}
                }
            }
        }
        if(length(x$names$competing)>0){
            for(k in competing_variables_position){
                later_variables=setdiff((k+1):NCOL(prepared_data),outcome_variables_position)
                later_outcome_variables=intersect((k+1):NCOL(prepared_data),outcome_variables_position)
                later_competing_variables=intersect((k+1):NCOL(prepared_data),competing_variables_position)
                if(any(has_died <- (prepared_data[[k]]%in%1))){
                    for(l in later_variables) {set(prepared_data,j=l,i=which(has_died),value=NA)}
                    for(l in later_outcome_variables) {set(prepared_data,j=l,i=which(has_died),value=0)}
                    for(l in later_competing_variables) {set(prepared_data,j=l,i=which(has_died),value=1)}
                }
            }
        }
        ## All variables should be NA as soon as censoring has occurred
        if(length(x$names$censoring)>0){
            for(k in censoring_variables_position){
                later_variables=(k+1):NCOL(prepared_data)
                if(any(has_censored <- (prepared_data[[k]]%in%x$names$censored_label))){
                    for(l in later_variables) {set(prepared_data,j=l,i=which(has_censored),value=NA)}
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
                later_variables=(k+1):NCOL(prepared_data)
                if(any(has_censored <- (prepared_data[[k]]%in%x$names$censored_label))){
                    for(l in later_variables) {set(prepared_data,j=l,i=which(has_censored),value=NA)}
                }
            }
        }
    }
    # Find last interval for individual where the individual is still 'at-risk' of experiencing
    # the outcome event and still 'uncensored'
    x <- count_followup(x = x,
                        prepared_data = prepared_data,
                        outcome_variables,
                        censoring_variables,
                        competing_variables,
                        max_time_horizon = max_time_horizon)
    x$prepared_data <- prepared_data[]
    x$names$name_time_covariates <- name_time_covariates
    x$names$name_baseline_covariates <- name_baseline_covariates
    x$names$name_constant_variables <- c(constant_variables,non_variables)
    x
}

######################################################################
### prepare_data.R ends here
