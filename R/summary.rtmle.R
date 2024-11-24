### summary.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 29 2024 (10:44) 
## Version: 
## Last-Updated: Nov 24 2024 (06:47) 
##           By: Thomas Alexander Gerds
##     Update #: 63
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' For each target in the object the risk estimates are summarized in a table with their
#' 95% confidence limits.
#'
#' If a target includes multiple protocols the results includes risk differences.
#' @title Summarizing the results of a register data analysis with the targeted minimum loss estimator 
#' @param object Object to be summarized
#' @param digits Number of decimals for the confidence intervals
#' @param targets Names of targets for which to compute the summary. Defaults to all targets in the object.
#' @param reference (Optional) Named list of reference protocols, one for each target. 
#' @param ... not used
#' @examples
#'
#' set.seed(112)
#' ld <- simulate_long_data(n = 181,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 2),register_format = TRUE)
#' x <- rtmle_init(intervals = 2,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
#' x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
#' add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
#' x <- long_to_wide(x,intervals = seq(0,2000,30.45*6))
#' protocol(x) <- list(name = "Always_A",treatment_variables = "A",intervention = 1)
#' protocol(x) <- list(name = "Never_A",treatment_variables = "A",intervention = 0)
#' prepare_data(x) <- list(reset = TRUE)
#' target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = c("Always_A","Never_A"))
#' x <- run_rtmle(x,time_horizon=1:2)
#' summary(x,digits=2)
#' 
#' @export
#' @method summary rtmle                 
summary.rtmle <- function(object,targets,reference = NULL,digits = 1,...){
    Estimate <- Upper <- Lower <- NULL
    if (missing(targets)) {
        targets <- names(object$targets)
    }
    else {
        stopifnot(all(targets %in% names(object$targets)))
    }
    sum <- do.call(rbind,lapply(targets,function(target_name){
        target <- object$targets[[target_name]]
        protocols <- target$protocols
        risk <- do.call(rbind,lapply(protocols,function(protocol_name){
            # to avoid the internal selfdetect problem we take a copy
            e <- data.table::copy(object$estimate[[target_name]][[protocol_name]])
            e[,"Estimate (CI_95)":= Publish::formatCI(x = 100*Estimate,lower = 100*Lower,upper = 100*Upper,show.x = TRUE,digits = digits)]
            e[]
        }))
        ## for (nix in c("Estimate","Standard_error","Lower","Upper"))
            ## set(risk,j = nix,value = NULL)
        # comparison
        if (length(protocols)>1){
            if (length(reference) == 0){
                ref <- protocols[[1]]
            }else{
                stopifnot(reference[[target_name]] %in% protocols)
                ref <- reference[[target_name]]
            }
            contrast <- do.call(rbind,lapply(setdiff(protocols,ref),function(protocol_name){
                do.call(rbind,lapply(unique(risk$Time_horizon),function(tp){
                    # FIXME: can the reference be taken out of the loop?
                    reference_estimate <- x$estimate[[target_name]][[ref]]$Estimate[[tp]]
                    reference_IC <- x$IC[[target_name]][[ref]][[tp]]
                    N <- NROW(reference_IC)
                    this_estimate <- x$estimate[[target_name]][[protocol_name]]$Estimate[[tp]]
                    this_IC <- x$IC[[target_name]][[protocol_name]][[tp]]
                    # risk difference
                    risk_difference_estimate = this_estimate - reference_estimate
                    risk_ratio_estimate = exp(log(this_estimate) - log(reference_estimate))
                    risk_difference_se = sd(this_IC - reference_IC)/sqrt(N)
                    risk_difference_lower = risk_difference_estimate - qnorm(.975)*risk_difference_se
                    risk_difference_upper = risk_difference_estimate + qnorm(.975)*risk_difference_se
                    risk_ratio_log_se <- sd(this_IC/this_estimate - reference_IC/reference_estimate)/sqrt(N)
                    risk_ratio_lower <- risk_ratio_estimate*exp(-qnorm(.975)*risk_ratio_log_se)
                    risk_ratio_upper <- risk_ratio_estimate*exp(qnorm(.975)*risk_ratio_log_se)
                    e <- data.table(Target = rep(target_name,2),
                               Protocol = rep(protocol_name, 2),
                               Target_parameter=c("ATE", "Risk ratio"),
                               Time_horizon = tp,
                               Estimator = x$estimate[[target_name]][[ref]]$Estimator,
                               Reference = rep(reference, 2),
                               Estimate = c(risk_difference_estimate, risk_ratio_estimate),
                               Standard_error = c(risk_difference_se, risk_ratio_log_se),
                               Lower = c(risk_difference_lower, risk_ratio_lower),
                               Upper = c(risk_difference_upper, risk_ratio_upper),
                               P_value = c(2*pnorm(-abs(risk_difference_estimate/risk_difference_se)), 2*pnorm(-abs(log(risk_ratio_estimate)/risk_ratio_log_se))))
                    e[Target_parameter == "ATE","Estimate (CI_95)":= Publish::formatCI(x = 100*Estimate,lower = 100*Lower,upper = 100*Upper,show.x = TRUE,digits = digits)]
                    e[Target_parameter == "Risk_Ratio","Estimate (CI_95)":= Publish::formatCI(x = Estimate,lower = Lower,upper = Upper,show.x = TRUE,digits = digits)]
                    e[]
                }))}))
            
            out <- data.table::rbindlist(list(risk,contrast),use.names = TRUE)
            return(out)
        }else{
            ## contrast <- NULL
            return(risk)
        }
    }))
    sum
}


######################################################################
### summary.rtmle.R ends here
