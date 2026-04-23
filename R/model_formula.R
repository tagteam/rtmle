### model_formula.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 16 2025 (08:58) 
## Version: 
## Last-Updated: apr 23 2026 (17:41) 
##           By: Thomas Alexander Gerds
##     Update #: 139
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
##' except at time 0. The reason for this is that we do not want to mistakenly assume that
##' L_k -> A_k when in reality A_k may happen before L_k.
##'
##' When interventions depend on more than one variable
##' @title Model formulas for nuisance parameters
##' @param x  object of class \code{rtmle}
##' @param propensity_model Only relevant for interventions which depend on multiple treatment variables:
##' Decide how to model the propensity of the regimen that the protocol dictates:
##' Possible values are \code{"joint"}, \code{"sequential"} or \code{"independent"}.
##' @param Markov Names of time-dependent variables which should only occur with the most recent
##'               value on the right hand side of the formulas.
##' @param exclude_variables Variables to exclude from the formulas for the nuisance parameters.
##' @param exclusion_rules Experimental. Additional exclusion rules given as a named list where names are variables that occur on the left hand side of a formula and
##'                        elements are variables that should be excluded from the right hand side of the formula. 
##' @param inclusion_rules Experimental. Additional inclusion rules given as a named list where names are variables that occur on the left hand side of a formula and
##' elements are variables that should be included in the right hand side of the formula. 
##' @param verbose Logical. If \code{FALSE} suppress all messages. \code{TRUE} is the default.
##' @param ... Not used (not yet)
##' @return The modified \code{rtmle}object
##' @examples
##' set.seed(17)
##' ld <- simulate_long_data(n = 17,number_visits = 20,
##'                          beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),
##'                          register_format = TRUE)
##' ld$timevar_data$B <- ld$timevar_data$A[id %in% sort(sample(1:91,size = 75,replace = FALSE))]
##' x <- rtmle_init(time_grid = seq(0,1500,30.45*12),name_id = "id",name_outcome = "Y",
##'                 name_competing = "Dead",
##'                 name_censoring = "Censored",censored_label = "censored")
##' x <- add_long_data(x,
##'                    outcome_data=ld$outcome_data,
##'                    censored_data=ld$censored_data,
##'                    competing_data=ld$competing_data,
##'                    timevar_data=ld$timevar_data)
##' x <- add_baseline_data(x,data=ld$baseline_data)
##' x <- long_to_wide(x,start_followup_date=0)
##' x <- prepare_rtmle_data(x)
##' x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1)
##' x <- protocol(x,name = "Never_A",treatment_variables = "A",intervention = 0)
##' x <- protocol(x,name = "Use_A_not_B",treatment_variables = c("A","B"),intervention = c(1,0))
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
                          propensity_model = "joint",
                          Markov = NULL,
                          verbose = TRUE,
                          exclude_variables = NULL,
                          exclusion_rules = NULL,
                          inclusion_rules = NULL,
                          ...){
    exclude_variables = c("start_followup_date",exclude_variables)
    name_constant_variables <- x$names$name_constant_variables
    if (length(x$protocols) == 0) {stop("No protocols registered in object, hence it is unclear which variables are intervened upon.")}
    all_treatment_variables <- setdiff(unlist(lapply(x$protocols,function(u)u$treatment_variables),use.names = FALSE),name_constant_variables)
    name_time_covariates <- setdiff(x$names$name_time_covariates,exclude_variables)
    name_baseline_covariates <- setdiff(x$names$name_baseline_covariates,exclude_variables)
    if (length(name_time_covariates)>0){
        if (length(Markov)>0 && Markov[[1]]!="")
            if (any(not_found <- !(Markov%in%name_time_covariates)))
                stop(paste0("The following variables in argument Markov do not match time_covariates:\n",
                            paste(Markov[not_found],collapse=", ")))
    }
    #
    # loop across time points looking at treatments, censoring, outcomes from the beginning of the interval
    # 
    model_formulas <- lapply(x$intervention_nodes,function(tk){
        all_vars <- c(lapply(x$protocols,function(pro){
            vals <- pro$intervention_table[time == tk][["value"]]
            if (length(vals) == 0){
                NULL
            }else{
                names(vals) <- pro$intervention_table[time == tk][["variable"]]
                vals
            }
        }))
        all_vars <- all_vars[!sapply(all_vars, is.null)]
        # censoring variables before outcome (no model for censoring at time zero)
        if(length(x$names$censoring)>0){
            censvalue <- paste0("'",x$names$uncensored_label,"'")
            names(censvalue) = paste0(x$names$censoring,"_",(tk+1))
            all_vars <- c(all_vars,list("censoring" = censvalue))
        }
        # outcome variables (no model for outcome at time zero)
        outvalue <- 1
        names(outvalue) <- paste0(x$names$outcome,"_",(tk+1))
        all_vars <- c(all_vars,list("outcome" = outvalue))
        # return vector of character formulas
        tk_forms <- lapply(names(all_vars), function(nav){
            vv <- all_vars[[nav]]
            ## this may be confusing but at intervention node k
            ## we evaluate treatment in the previous interval [t_{k-1},t_k] but
            ## censoring and outcome in the next interval [t_{k},t_{k+1}]
            if (nav%in% c("censoring","outcome"))
                eval_time <- tk+1
            else
                eval_time <- tk
            ff <- formalize(timepoint = eval_time,
                            available_names = names(x$prepared_data),
                            name_outcome_variable = names(vv),
                            outcome_value = as.character(vv),
                            name_baseline_covariates = name_baseline_covariates,
                            name_time_covariates  = name_time_covariates,
                            Markov = Markov,
                            constant_variables = name_constant_variables,
                            exclusion_rules = exclusion_rules,
                            inclusion_rules = inclusion_rules,
                            handle_concomitant_variables = propensity_model,
                            unwanted_variables = exclude_variables)
            ff
        })
        names(tk_forms) <- names(all_vars)
        tk_forms
    })
    names(model_formulas) <- paste0("time_",x$intervention_nodes)
    x$models <- model_formulas
    x
}


######################################################################
### model_formula.R ends here
