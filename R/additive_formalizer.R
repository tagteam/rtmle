### additive_formalizer.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (11:13)
## Version:
## Last-Updated: Jul 24 2025 (09:44) 
##           By: Thomas Alexander Gerds
##     Update #: 96
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
                                Markov = NULL) {
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
    # loop across time points
    propensity_formulas <- c(unlist(lapply(x$times,function(tk){
        c(unlist(lapply(paste0(treatment_variables,"_",tk), function(reg){
            formalize(timepoint = tk,
                      available_names = names(x$prepared_data),
                      name_outcome_variable = reg,
                      name_baseline_covariates = name_baseline_covariates,
                      name_time_covariates  = name_time_covariates,
                      Markov = Markov,
                      constant_variables = name_constant_variables)
        })))
    })))
    if(length(x$names$censoring)>0){
        censoring_formulas <- c(unlist(lapply(x$times[-1],function(tk){
            formalize(timepoint = tk,
                      available_names = names(x$prepared_data),
                      name_outcome_variable = paste0(x$names$censoring,"_",tk),
                      name_baseline_covariates = name_baseline_covariates,
                      name_time_covariates = name_time_covariates,
                      Markov = Markov,
                      constant_variables = name_constant_variables)
        })))
    }
    ## Note that A_k ~ V + L_0 + ... + L_(k-1) + A_(k-1) for k = 1,..., max(x$times), but A_0 ~ V + L_0
    ## i.e., regimen at baseline depends on additional baseline covariates, whereas in general, regimen depends
    ## on the previously observed covariates and regimen.
    ## The reason for this is we do not want to mistakenly assume that L_1 -> A_1 when in reality A_1 happens before L_1
    outcome_formulas <- unlist(lapply(x$times[-1],function(tk){
        formalize(timepoint = tk,
                  available_names = names(x$prepared_data),
                  name_outcome_variable = paste0(x$names$outcome,"_",tk),
                  name_baseline_covariates = name_baseline_covariates,
                  name_time_covariates  = name_time_covariates,
                  Markov = Markov, constant_variables = name_constant_variables)
    }))
    names(outcome_formulas)=paste0(x$names$outcome,"_",x$times[-1])
    # names for treatment and censoring formulas
    names(propensity_formulas) <- as.character(unlist(lapply(propensity_formulas,function(x){strsplit(x," ~ ")[[1]][[1]]})))
    if (length(x$names$censoring)>0){
        names(censoring_formulas) <- as.character(unlist(lapply(censoring_formulas,function(x){strsplit(x," ~ ")[[1]][[1]]})))
        # Remove gformulas of variables not in data, e.g., if we exclude censor others at time point 0
        censoring_formulas <- censoring_formulas[names(censoring_formulas) %in% names(x$prepared_data)]
    }else{
        censoring_formulas <- NULL
    }
    # Remove gformulas of variables not in data, e.g., if we exclude censor others at time point 0
    propensity_formulas <- propensity_formulas[names(propensity_formulas) %in% names(x$prepared_data)]
    lapply(c(propensity_formulas,censoring_formulas,outcome_formulas),function(f)list(formula = f))
}
######################################################################
### additive_formalizer.R ends here
