### formalize.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  4 2024 (07:40) 
## Version: 
## Last-Updated: mar 16 2026 (15:35) 
##           By: Thomas Alexander Gerds
##     Update #: 90
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
                      outcome_value,
                      name_baseline_covariates,
                      name_time_covariates,
                      Markov = NULL,
                      constant_variables = NULL,
                      exclusion_rules = NULL,
                      inclusion_rules = NULL,
                      handle_concomitant_variables,
                      unwanted_variables = NULL){
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
    # remove outcome variable (could be removed earlier in the program)
    included_vars <- setdiff(included_vars,name_outcome_variable)
    # remove unwanted variables
    included_vars <- setdiff(included_vars,unwanted_variables)    
    # remove vars that are not in data (if data are available yet)
    included_vars <- intersect(included_vars,available_names)
    # apply exclusion_rules
    if (length(exclusion_rules)>0){
        matches <- vapply(names(exclusion_rules), function(e) {grepl(e, paste0(name_outcome_variable,collapse = " "))}, logical(1))
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
    if (length(inclusion_rules)>0){
        matches <- vapply(names(inclusion_rules), function(e) {grepl(e, paste0(name_outcome_variable,collapse = " "))}, logical(1))
        if (sum(matches) > 0){
            inclusion_pattern <- c(unlist(inclusion_rules[matches]))
            included_additions <- grep(paste0(inclusion_pattern,collapse = "|"),available_names,value = TRUE)
            if (length(included_additions)>0){
                included_vars <- unique(c(included_vars,inclusion_pattern))
            }else{
                warning(paste0("The following inclusion_rule did not match any variable in the data:\n",paste0(names(matches),collapse = ",")))
            }
        }
    }
    # In any case the outcome variable(s) cannot be predictor(s)
    included_vars <- setdiff(included_vars,name_outcome_variable)
    if(length(included_vars)>0){
        rhs = paste(included_vars, collapse = " + ")
    }else{
        rhs <- "1"
    }
    if ((nvars <- length(name_outcome_variable))>1){
        if (handle_concomitant_variables == "joint"){
            outcome_string <- paste0("I(1*(",paste0(paste0(name_outcome_variable,"==",outcome_value),collapse = "&"),"))")
        }else{
            outcome_string <- paste0("I(1*(",name_outcome_variable,"==",outcome_value,"))")
        }
        if (handle_concomitant_variables == "sequential"){
            form <- vector(mode = "list",length = nvars)
            for (j in 1:nvars){
                concomitant_rhs <- paste0(name_outcome_variable[-c(j:nvars)],collapse = "+")
                if (rhs != "1"){
                    concomitant_rhs <- paste0(rhs,"+",concomitant_rhs)
                }else{
                    concomitant_rhs <- "1"
                }
                form[[j]] <- list(list(formula = paste0(outcome_string[[j]]," ~ ",concomitant_rhs)))
                names(form[[j]]) <- name_outcome_variable[[j]]
            }
        }else{
            # joint 
            if (handle_concomitant_variables == "joint"){
                form <- list(list(formula = paste0(outcome_string," ~ ",rhs)))
                names(form) <- paste0(name_outcome_variable,collapse = ",")
            }else{
                # independent
                form <- do.call("c",lapply(outcome_string,function(os){
                    list(list(formula = paste0(os," ~ ",rhs)))
                }))
                names(form) <- name_outcome_variable
            }
        }
    }else{
        outcome_string <- paste0("I(1*(",name_outcome_variable,"==",outcome_value,"))")
        form <- list(list(formula = paste0(outcome_string," ~ ",rhs)))
        names(form) <- name_outcome_variable
    }
    form
}
######################################################################
### formalize.R ends here
