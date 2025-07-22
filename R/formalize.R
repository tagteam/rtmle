### formalize.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  4 2024 (07:40) 
## Version: 
## Last-Updated: Jul 22 2025 (09:20) 
##           By: Thomas Alexander Gerds
##     Update #: 26
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
formalize <- function(timepoint,
                      work_data,
                      name_outcome_variable,
                      name_baseline_covariates,
                      name_time_covariates,
                      Markov = NULL,
                      constant_variables = NULL,
                      exclusion_rules = NULL){
    form = NULL
    # remove constant variables
    included_baseline_covariates <- setdiff(name_baseline_covariates,constant_variables)
    # check Markov assumption for time-varying variables
    # and add _time 
    has_markov=match(name_time_covariates,Markov,nomatch=0)
    if (sum(has_markov)>0){
        markov_time_covariates=sapply(name_time_covariates[has_markov>0],function(v){paste0(v,"_",max(0,timepoint-1))})
    } else{
        markov_time_covariates=NULL
    }
    if (any(has_markov==0)){
        non_markov_time_covariates=sapply(name_time_covariates[has_markov==0],function(v){paste0(v,"_",0:max(0,timepoint-1))})
    }    else{
        non_markov_time_covariates=NULL
    }
    included_vars=c(included_baseline_covariates,
                    setdiff(non_markov_time_covariates,constant_variables),
                    setdiff(markov_time_covariates,constant_variables))
    # remove outcome variable (could be removed in earlier functions)
    included_vars = setdiff(included_vars,name_outcome_variable)
    # remove vars that are not in data
    has_not=match(included_vars,names(work_data),nomatch=0)==0
    included_vars=included_vars[!has_not]
    # apply exclusion_rules
    if (length(exclusion_rules)>0 && name_outcome_variable %in% names(exclusion_rules))
    included_vars <- setdiff(included_vars,exclusion_rules[[name_outcome_variable]])
    # add covariates if any
    if(length(included_vars)>0){
        form = paste(included_vars, collapse = " + ")
    }else{
        form <- "1"
    }
    # return character formula
    form = paste0(name_outcome_variable," ~ ",form)
    form[]
}
######################################################################
### formalize.R ends here
