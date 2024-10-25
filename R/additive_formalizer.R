### additive_formalizer.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (11:13) 
## Version: 
## Last-Updated: Oct 22 2024 (09:44) 
##           By: Thomas Alexander Gerds
##     Update #: 65
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
additive_formalizer <- function(x,
                                protocol,
                                treatment_variables,
                                include_variables = NULL,
                                exclude_variables = NULL,
                                Markov = NULL){
    # FIXME: improve the stopping message
    if (missing(treatment_variables))
        treatment_variables <- x$protocols[[protocol]]$treatment_variables
    stopifnot(length(treatment_variables)>0)
    name_time_covariates <- setdiff(x$names$name_time_covariates,exclude_variables)
    name_baseline_covariates <- setdiff(x$names$name_baseline_covariates,exclude_variables)
    name_constant_variables <- x$names$name_constant_variables
    if (length(name_time_covariates)>0){
        if (length(Markov)>0 && Markov[[1]]!="")
            if (any(not_found <- !(Markov%in%name_time_covariates)))
                stop(paste0("The following variables in argument Markov do not match time_covariates:\n",
                            paste(Markov[not_found],collapse=", ")))
    }
    K = length(x$times)
    # baseline propensity
    propensity_formulas <- c(unlist(lapply(treatment_variables,function(reg){
        paste0(reg,"_0"," ~ ",formalize(timepoint = 0,
                                        work_data = x$prepared_data,
                                        name_baseline_covariates = name_baseline_covariates,
                                        name_time_covariates = name_time_covariates,
                                        Markov = Markov,
                                        # remove baseline treatment from predictor variables
                                        constant_variables = c(name_constant_variables,reg)))
    })))
    # censoring in first interval
    if(length(x$names$censoring)>0){
        censoring_formulas <- paste0(x$names$censoring,"_1"," ~ ", formalize(timepoint = 0,
                                                                             work_data = x$prepared_data,
                                                                             name_baseline_covariates = name_baseline_covariates,
                                                                             name_time_covariates = name_time_covariates,
                                                                             Markov = Markov,
                                                                             constant_variables = name_constant_variables))
    } else {
        censoring_formulas <- NULL
    }
    # other intervals
    if(K>1){
        propensity_formulas <- c(propensity_formulas,unlist(lapply(x$times[-1],function(tk){
            c(unlist(lapply(treatment_variables, function(reg){
                paste0(reg,"_",tk," ~ ",formalize(timepoint = tk,
                                                  work_data = x$prepared_data,
                                                  name_baseline_covariates = name_baseline_covariates,
                                                  name_time_covariates  = name_time_covariates,
                                                  Markov = Markov,
                                                  constant_variables = name_constant_variables))
            })))
        })))
        if(length(x$names$censoring)>0){
            censoring_formulas <- c(censoring_formulas,unlist(lapply(x$times[-c(1,2)],function(tk){
                paste0(x$names$censoring,"_",tk," ~ ",formalize(timepoint = tk,
                                                                work_data = x$prepared_data,
                                                                name_baseline_covariates = name_baseline_covariates,
                                                                name_time_covariates = name_time_covariates,
                                                                Markov = Markov,
                                                                constant_variables = name_constant_variables))
            })))
        }
    }
    ## Note that A_k ~ V + L_0 + ... + L_(k-1) + A_(k-1) for k = 1,..., max(x$times), but A_0 ~ V + L_0
    ## i.e., regimen at baseline depends on additional baseline covariates, whereas in general, regimen depends
    ## on the previously observed covariates and regimen.
    ## The reason for this is we do not want to mistakenly assume that L_1 -> A_1 when in reality A_1 happens before L_1
    outcome_formulas <- unlist(lapply(x$times[-1],function(tk){
        paste0(x$names$outcome,"_",tk," ~ ", formalize(timepoint = tk, work_data = x$prepared_data,
                                                       name_baseline_covariates = name_baseline_covariates,
                                                       name_time_covariates  = name_time_covariates, 
                                                       Markov = Markov, constant_variables = name_constant_variables))
    }))
    names(outcome_formulas)=paste0(x$names$outcome,"_",x$times[-1])
    # names for treatment and censoring formulas
    names(propensity_formulas) <- as.character(unlist(lapply(propensity_formulas,function(x){strsplit(x," ~ ")[[1]][[1]]})))
    names(censoring_formulas) <- as.character(unlist(lapply(censoring_formulas,function(x){strsplit(x," ~ ")[[1]][[1]]})))
    # Remove gformulas of variables not in data, e.g., if we exclude censor others at time point 0
    propensity_formulas <- propensity_formulas[names(propensity_formulas) %in% names(x$prepared_data)]
    censoring_formulas <- censoring_formulas[names(censoring_formulas) %in% names(x$prepared_data)]
    lapply(c(propensity_formulas,censoring_formulas,outcome_formulas),function(f)list(formula = f))
}
######################################################################
### additive_formalizer.R ends here
