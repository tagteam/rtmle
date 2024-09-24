### formalize.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  4 2024 (07:40) 
## Version: 
## Last-Updated: Sep 23 2024 (12:11) 
##           By: Thomas Alexander Gerds
##     Update #: 15
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
                      name_baseline_covariates,
                      name_time_covariates,
                      Markov = NULL,
                      constant_variables = NULL){
    form = NULL
    # remove constant variables
    included_baseline_covariates=setdiff(name_baseline_covariates,constant_variables)
    included_time_covariates=setdiff(name_time_covariates,constant_variables)
    # check Markov assumption for time-varying variables
    # and add _time 
    has_markov=match(included_time_covariates,Markov,nomatch=0)
    if (sum(has_markov)>0){
        markov_time_covariates=sapply(included_time_covariates[has_markov>0],function(v){paste0(v,"_",max(0,timepoint-1))})
    } else{
        markov_time_covariates=NULL
    }
    if (any(has_markov==0)){
        non_markov_time_covariates=sapply(included_time_covariates[has_markov==0],function(v){paste0(v,"_",0:max(0,timepoint-1))})
    }    else{
        non_markov_time_covariates=NULL
    }
    included_vars=c(included_baseline_covariates,
                    non_markov_time_covariates,
                    markov_time_covariates)
    # remove vars that are not in data
    has_not=match(included_vars,names(work_data),nomatch=0)==0
    included_vars=included_vars[!has_not]
    # add covariates if any
    if(length(included_vars)>0){
        form = paste(included_vars, collapse = " + ")
    }else{
        form <- "1"
    }
    form[]
}
######################################################################
### formalize.R ends here
