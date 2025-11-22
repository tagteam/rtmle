### model_formula.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 16 2025 (08:58) 
## Version: 
## Last-Updated: nov 20 2025 (09:17) 
##           By: Thomas Alexander Gerds
##     Update #: 58
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Specify model formulas for the estimators of the nuisance parameters
##' 
##'
##' Note that due to the discretized nature of the data, we do not include
##' the time-dependent covariates at time k in the models for treatment at time k,
##' except at time 0. 
##' The reason for this is we do not want to mistakenly assume that L_k -> A_k when in reality A_k may happen before L_k
##' @title Model formulas for nuisance parameters
##' @param x  object of class \code{rtmle}
##' @param Markov Names of time-dependent variables which should only occur with the most recent
##'               value on the right hand side of the formulas.
##' @param exclude_variables Variables to exclude from the formulas for the nuisance parameters.
##' @param exclusion_rules Experimental. Additional exclusion rules given as a named list where names are variables that occur on the left hand side of a formula and
##'                        elements are variables that should be included in the right hand side of the formula. 
##' @param inclusion_rules Experimental. Additional inclusion rules given as a named list where names are variables that occur on the left hand side of a formula and
##' elements are variables that should be included in the right hand side of the formula. 
##' @param verbose Logical. If \code{FALSE} suppress all messages. \code{TRUE} is the default.
##' @param ... Not used (not yet)
##' @return The modified \code{rtmle}object
##' @examples
##' set.seed(112)
##' ld <- simulate_long_data(n = 1000,number_visits = 20,
##'                          beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),
##'                          register_format = TRUE)
##' x <- rtmle_init(intervals = 2,name_id = "id",name_outcome = "Y",
##'                 name_competing = "Dead",
##'                 name_censoring = "Censored",censored_label = "censored")
##' x <- add_long_data(x,outcome_data = ld[["outcome_data"]],
##'                      censored_data = ld[["censored_data"]],
##'                      competing_data = ld[["competing_data"]],
##'                      timevar_data = ld[["timevar_data"]])
##' x <- add_baseline_data(x,data = ld$baseline_data)
##' x <- long_to_wide(x,intervals = seq(0,2000,30.45*6))
##' x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1)
##' x <- prepare_data(x) 
##' x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "Always_A")
##' x <- model_formula(x)
##' x$models
##' # remove age from all formulas
##' x <- model_formula(x,exclusion_rules=list("*"="age"))
##' # remove age from all Y_1 formula
##' x <- model_formula(x,exclusion_rules=list("Y_1"="age"))
##' x$models
##' # remove age from all Y formula
##' x <- model_formula(x,exclusion_rules=list("Y_*"="age"))
##' # remove age and L_0 from Censored_1 and from Censored_2 formula
##' x <- model_formula(x,exclusion_rules=list("Censored_[1:2]"="age|L_0"))
##' x$models
##' # remove all L_t variables from Censored_1 and from Censored_2 formula
##' x <- model_formula(x,exclusion_rules=list("Censored_[1:2]"="L_*"))
##' x$models
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
model_formula <- function(x,
                          Markov = NULL,
                          verbose = TRUE,
                          exclude_variables = NULL,
                          exclusion_rules = NULL,
                          inclusion_rules = NULL,
                          ...){
    exclude_variables = c("start_followup_date",exclude_variables)
    name_constant_variables <- x$names$name_constant_variables
    all_treatment_variables <- setdiff(unlist(lapply(x$protocols,function(u)u$treatment_variables),use.names = FALSE),name_constant_variables)
    name_time_covariates <- setdiff(x$names$name_time_covariates,exclude_variables)
    name_baseline_covariates <- setdiff(x$names$name_baseline_covariates,exclude_variables)
    if (length(name_time_covariates)>0){
        if (length(Markov)>0 && Markov[[1]]!="")
            if (any(not_found <- !(Markov%in%name_time_covariates)))
                stop(paste0("The following variables in argument Markov do not match time_covariates:\n",
                            paste(Markov[not_found],collapse=", ")))
    }
    # loop across time points
    model_formulas <- c(unlist(lapply(x$times,function(tk){
        if (tk >= (length(x$times)-1)){ # no treatment formula for last time point, minus 1 because 0 is included in x$times
            all_vars <- NULL
        }else{
            all_vars <- unlist(lapply(all_treatment_variables,function(tv)paste0(tv,"_",tk)))
        }
        if (tk != 0){
            # outcome variables (no model for outcome at time zero)
            all_vars <- c(all_vars,paste0(x$names$outcome,"_",tk))
            # censoring variables (no model for censoring at time zero)
            if(length(x$names$censoring)>0){    
                all_vars <- c(all_vars,paste0(x$names$censoring,"_",tk))
            }
        }
        # remove constant variables 
        # FIXME: this also removes constant outcome variables
        ## all_vars <- setdiff(all_vars,name_constant_variables)
        # return vector of character formulas
        c(unlist(lapply(all_vars, function(vv){
            ff <- formalize(timepoint = tk,
                            available_names = names(x$prepared_data),
                            name_outcome_variable = vv,
                            name_baseline_covariates = name_baseline_covariates,
                            name_time_covariates  = name_time_covariates,
                            Markov = Markov,
                            constant_variables = name_constant_variables,
                            exclusion_rules = exclusion_rules,
                            inclusion_rules = inclusion_rules,
                            unwanted_variables = exclude_variables)
            names(ff) = vv
            ff
        })))
    })))
    # Remove formulas of variables that do not occur in data
    if (verbose && any(do_not_occur <- !(names(model_formulas) %in% names(x$prepared_data)))){
        message(paste0("The right hand sides of the following formulas do not occur in the data, hence we drop them:\n", 
                       paste0(model_formulas[do_not_occur],collapse = "\n")))
    }
    model_formulas <- model_formulas[names(model_formulas) %in% names(x$prepared_data)]
    x$models <- lapply(model_formulas,function(f)list(formula = f))
    if (verbose)
    x
}


######################################################################
### model_formula.R ends here
