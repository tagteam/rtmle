### obsolete_functions.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 16 2025 (14:50) 
## Version: 
## Last-Updated: dec  3 2025 (06:21) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
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
"prepare_data<-" <- function(x,...,value){
    .Deprecated("prepare_data")
    # check value
    stopifnot(is.list(value))

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
    check_id <- function(data,id,allow_null = TRUE, name = ""){
        if (NROW(data)>0){
            if (!inherits(data,"data.table"))
                if (inherits(data,"data.frame")){
                    data.table::setDT(data)
                } else{
                    warning("The object x$outcome_data must be a data.table, data.frame, or tibble")
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
    nix = lapply(names(x$data$timevar_data),function(nn){
        check_id(data = x$data$timevar_data[[nn]], id = x$names$id,name = paste0("x$data$timevar_data$",nn))
    })
    rm(nix)
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
        name_baseline_covariates = setdiff(names(x$data$baseline_data),x$names$id)
        if (!is.null(x$data$baseline_data)){
            stopifnot(inherits(x$data$baseline_data,"data.table"))
            stopifnot(match(x$names$id,names(x$data$baseline_data),nomatch = 0)>0)
            prepared_data=x$data$baseline_data[prepared_data,on = x$names$id]
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
            prepared_data <- prepared_data[x$data$timevar_data,on = x$names$id]
        }else{
            name_time_covariates <- names(x$data$timevar_data)
            if (length(name_time_covariates)>0){
                for (Vname in name_time_covariates){
                    prepared_data <- prepared_data[x$data$timevar_data[[Vname]],on = x$names$id]
                }
            }
        }
    } else {
        name_time_covariates <- NULL
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
            } else
                paste0(c(x$names$censoring, x$names$outcome),"_",timepoint)
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
    x$names$name_constant_variables <- constant_variables
    x
}

#' Adding baseline data to a rtmle object
#' 
#' This function adds a baseline dataset to an existing rtmle object
#' @param x object of class \code{rtmle} 
#' @param ... Not used (not yet)
#' @param value data.frame, data.table or tibble which contains the baseline data
#'@export
"add_baseline_data<-" <- function(x,...,value){
    .Deprecated("add_baseline_data")
    if (!is.data.frame(value) || !(x$name$id %in% names(value)))
        stop(paste0("'value' needs to be a data.frame (or data.table or tibble) with a variable called ",x$names$id))
    if (!(x$names$id %in% names(value)))
        stop(paste0("'value' needs to contain the subject identifier variable: '",x$names$id,"'"))
    data.table::setDT(value)
    x$data$baseline_data <- value
    x
}


##' Define a named protocol for a hypothetical/emulated trial
##'
##' This function adds a protocol to an existing object.
##' A protocol defines the values of the treatment variable(s)
##' at each time point during followup including at time zero (baseline).
##' @param x An rtmle object as obtained with \code{rtmle_init}.
##' @param ... Not (yet) used
##' @param value list with three forced elements:
##' \itemize{
##' \item \code{name}: the name of the protocol
##' \item \code{treatment_variables}: A vector with the
##' name(s) of the variable(s) that the protocols intervenes upon. In longitudinal 
##' settings, when a treatment variable is named "A" then the prepared data will contain
##' a column for each of the variables A_0,A_1,----,A_k. This argument can be left unspecified
##' in which case the argument \code{intervention} must be a data.frame with names of treatment variables.
##' An optional argument is
##' \itemize{
##' \item \code{verbose}: if FALSE suppress messages
##' }
##' See also details.
##' \item \code{intervention}: A vector, or a named data.frame (tibble, data.table), or a function.
##' If it is a vector, it should consist of values 0 and 1 corresponding to what the intervention
##' would set for the \code{treatment_variables}. If it is a data.frame (tibble, data.table),
##' the names specify the treatment variables and the argument \code{treatment_variables} is ignored.
##' Each column must be a factor with levels specifying the treatment options. In longitudinal settings,
##' there can either be only one row in the data.frame and then it is assumed that the intervention is
##' static throughout the followup period. The data.frame can also have a column named time and
##' specify for each time point the values that the intervention sets for the treatment variables.
##' If it is a function, it will be called from \code{intervention_probabilities} with two arguments:
##' the current time interval and the current history of all variables. The function determines the value(s)
##' of the treatment variable(s) under the intervention and should return a matrix with as many columns as there are
##' treatment variables.
##' }
#' @return The modified object contains the treatment variables and the intervention_table
#'         as list elemens of \code{x$protocols[[name]]} where name is given by \code{value$name}
#' @author  Thomas A Gerds \email{tag@@biostat.ku.dk} 
#' @examples
#' # ------------------------------------------------------------------------------------------
#' # Intervening on a single treatment variable
#' # ------------------------------------------------------------------------------------------
#' x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- protocol(x,name = "Always_A",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1"))))
#' x <- protocol(x,name = "Never_A",
#'                     intervention = data.frame("A" = factor("0",levels = c("0","1"))))
#' x <- protocol(x,name = "Initiate_A_then_stop",
#'                     intervention = data.frame("A" = factor(c("1","0","0"),levels = c("0","1"))))
#' x$protocols
#' # ------------------------------------------------------------------------------------------
#' # Intervening on a more than one treatment variable
#' # ------------------------------------------------------------------------------------------
#' x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- protocol(x,name = "Always_A_never_B",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1")),
#'                                               "B" = factor("0",levels = c("0","1"))))
#' x <- protocol(x,name = "Always_A_and_B_never_C",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1")),
#'                                               "B" = factor("1",levels = c("0","1")),
#'                                               "C" = factor("0",levels = c("0","1"))))
#' x$protocols
##' @export
"protocol<-" <- function(x,...,value) {
    .Deprecated("protocol")
    variable <- NULL
    stopifnot(is.list(value))
    stopifnot(all(c("name","intervention") %in% names(value)))
    if (length(value$verbose) == 0) value$verbose <- TRUE
    intervention_times <- x$time[-length(x$time)]
    intervention_table <- value$intervention
    if (inherits(intervention_table,"data.frame")){
        data.table::setDT(intervention_table)
        treatment_variables <- names(intervention_table)
        if (any(grepl(pattern = "_[0-9]+$",
                      x = treatment_variables))){
            stop("Treatment variables should be given without time suffix.")
        }
        missing_nodes <- (length(intervention_times)-NROW(intervention_table))
        if (missing_nodes>0){
            if (value$verbose[[1]] == TRUE){
                message("The object specifies more intervention nodes than there are rows in the provided intervention table.\nApply last value carried forward for now, 'x$protocol$intervention_table'.")
            }
            intervention_table <- intervention_table[c(1:NROW(intervention_table),rep(NROW(intervention_table),missing_nodes))]
        }else{
            if (missing_nodes<0){
                if (value$verbose[[1]] == TRUE){
                    message("The object specifies fewer intervention nodes than there are rows in the provided intervention table.\nCutting these for now, but please check 'x$protocol$intervention_table'.")
                }
                intervention_table <- intervention_table[c(1:length(intervention_times))]
            }
        }
        treatment_options <- sapply(treatment_variables,
                                    function(v){
                                        if (is.factor(value$intervention[[v]])){
                                            levels(value$intervention[[v]])
                                        } else{
                                            # FIXME: this seems out of control and will not work
                                            stop("The treatment variables must be factors") 
                                            unique(value$intervention[[v]])
                                        }
                                    },simplify = FALSE)
    }else{
        if ("treatment_variables" %in% names(value) &&
            length(value$treatment_variables) == length(intervention_table) &&
            all(intervention_table %in% c(0,1))){
            treatment_variables <- value$treatment_variables
            intervention_table <- data.table::as.data.table(lapply(1:length(treatment_variables),function(v){
                factor(intervention_table[[v]],levels = c(0,1))
            }))
            data.table::setnames(intervention_table,treatment_variables)
            treatment_options <- sapply(treatment_variables,function(x)c(0,1),simplify = FALSE)
        }else{
            stop("Intervention must be a vector of 0s and 1s or a data.frame (or data.table or tibble) containing factors.")
        }
    }
    if (!("time" %in% names(intervention_table))){
        intervention_table <- cbind(time = intervention_times,
                                    intervention_table)
    }else{
        if (!(all(intervention_table[["time"]] == intervention_times))){
            stop(paste("Intervention times do not match. Object contains:\n",
                       paste(intervention_times,collapse = ", ")))
        }
    }
    intervention_table <- do.call("rbind",lapply(treatment_variables,function(v){
        itv <- intervention_table[,c("time",v),with = FALSE]
        itv[,variable := paste0(v,"_",itv$time)]
        data.table::setnames(itv,old = v,new = "value")
        data.table::setcolorder(itv,c("time","variable","value"))
        itv
    }))
    if (length(value$intervene_function)>0){
        x$protocols[[value$name]]$intervene_function <- value$intervene_function
    }else{
        x$protocols[[value$name]]$intervene_function <- "intervene"
    }
    x$protocols[[value$name]]$treatment_variables <- treatment_variables
    x$protocols[[value$name]]$intervention_table <- intervention_table
    # adding the treatment options if necessary
    if (length(x$names$treatment_options) == 0){
        x$names$treatment_options <- treatment_options
    }else{
        new_options <- setdiff(names(treatment_options),
                               names(x$names$treatment_options))
        if (length(new_options)>0){
            x$names$treatment_options <- c(x$names$treatment_options,treatment_options[new_options])
        }
    }
    x
}

##' Define a target parameter 
##'
##' A target defines a target parameter, an estimator and models
##' for the nuisance parameters.
##' @param x An rtmle object as obtained with \code{rtmle_init}.
##' @param ... Not (yet) used
##' @param value list with three forced elements:
##' \itemize{
##' \item \code{name}: the name of the target parameter.
##' \item \code{strategy}: the nuisance parameter modeling strategy consists of formulas and libraries for the regression models.
##' \item \code{estimator}: the estimator of the target parameter.
##' }
##' @export
"target<-" <- function(x,...,value) {
    .Deprecated("target")
    stopifnot(is.list(value))
    stopifnot(all(c("name","strategy","estimator","protocols")%in%names(value)))
    stopifnot(all(value[["protocols"]]%in%names(x$protocols)))
    if (length(value[["include_variables"]])>0)
        include_variables <- value[["include_variables"]]
    else
        include_variables <- NULL
    if (length(value[["exclude_variables"]])>0)
        exclude_variables <- value[["exclude_variables"]]
    else
        exclude_variables <- NULL
    all_treatment_variables <- c(sapply(x$protocols,function(u)u$treatment_variables))
    if (length(value$strategy) == 1 && value$strategy == "additive"){
        # FIXME: some formulas could be shared across protocols
        for (protocol in value$protocols){
            if (length(value[["exclude_other_treatments"]])>0){
                protocol_exclude_variables <- c(exclude_variables,setdiff(all_treatment_variables,x$protocols[[protocol]]$treatment_variables))
            } else{
                protocol_exclude_variables <- exclude_variables
            }
            x$models[[protocol]] = additive_formalizer(x = x,
                                                       protocol = protocol,
                                                       exclude_variables = protocol_exclude_variables,
                                                       include_variables = include_variables,
                                                       Markov = value$markov)
            ## model(x) <- list(formalizer = "additive",treatment_variables = x$protocols[[value$protocol]]$treatment_variables)
            x$targets[[value$name]][["strategy"]] <- "additive"
        }
    }else{
        stop("Don't know about this strategy.")
    }
    x$targets[[value$name]][["protocols"]] <- unique(c(x$targets[[value$name]][["protocols"]],value$protocols))
    x$targets[[value$name]][["estimator"]] <- value$estimator
    x
}


######################################################################
### obsolete_functions.R ends here
