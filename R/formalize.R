### formalize.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  4 2024 (07:40) 
## Version: 
## Last-Updated: nov 20 2025 (09:15) 
##           By: Thomas Alexander Gerds
##     Update #: 54
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
formalize <- function(timepoint,
                      available_names,
                      name_outcome_variable,
                      name_baseline_covariates,
                      name_time_covariates,
                      Markov = NULL,
                      constant_variables = NULL,
                      exclusion_rules = NULL,
                      inclusion_rules = NULL,
                      unwanted_variables = NULL){
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
        non_markov_time_covariates <- sapply(name_time_covariates[has_markov==0],function(v){paste0(v,"_",0:max(0,timepoint-1))})
    }    else{
        non_markov_time_covariates <- NULL
    }
    included_vars <- c(included_baseline_covariates,
                       setdiff(non_markov_time_covariates,constant_variables),
                       setdiff(markov_time_covariates,constant_variables))
    # remove outcome variable (could be removed in earlier functions)
    included_vars <- setdiff(included_vars,name_outcome_variable)
    # remove unwanted variables
    included_vars <- setdiff(included_vars,unwanted_variables)    
    # remove vars that are not in data
    included_vars <- intersect(included_vars,available_names)
    # apply exclusion_rules
    if (length(exclusion_rules)>0){
        matches <- vapply(names(exclusion_rules), function(e) {grepl(e, name_outcome_variable)}, logical(1))
        if (sum(matches) > 0){
            exclusion_pattern <- c(unlist(exclusion_rules[matches]))
            # Regex metacharacters . \ | ( ) [ { ^ $ * + ? 
            ## regexp_pattern <- "[\\.\\\\\\|\\(\\)\\[\\{\\^\\$\\*\\+\\?]"
            ## has_regex <- grepl(regexp_pattern, exclusion_pattern)
            excluded_vars <- grep(paste0(exclusion_pattern,collapse = "|"),included_vars,value = TRUE)
            included_vars <- setdiff(included_vars,excluded_vars)
        }
    }
    # apply inclusion_rules
    if (length(inclusion_rules)>0 && name_outcome_variable %in% names(inclusion_rules)){
        if (any(has_not <- (match(inclusion_rules[[name_outcome_variable]],available_names,nomatch = 0) == 0))){
            warning(paste0("The following variables from the inclusion_rules do not occur in the data:\n",
                           paste0(inclusion_rules[[name_outcome_variable]][has_not],collapse = ", ")))
        }
        ivars <- intersect(inclusion_rules[[name_outcome_variable]],available_names)
        if (length(ivars)>0) included_vars <- unique(c(included_vars,ivars))
    }
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
