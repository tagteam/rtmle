### summary.rtmle.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 29 2024 (10:44) 
## Version: 
## Last-Updated: May 22 2025 (12:45) 
##           By: Thomas Alexander Gerds
##     Update #: 175
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
#' @param analysis Name of the analysis to be summarized
#' @param digits Number of decimals for the confidence intervals
#' @param targets Names of targets for which to compute the summary. Defaults to all targets in the object.
#' @param reference (Optional) Named list of reference protocols, one for each target.
#' @param ... not used
#' @examples
#'
#' set.seed(112)
#' ld <- simulate_long_data(n = 181,number_visits = 20,
#'                          beta = list(A_on_Y = -.2,
#'                                      A0_on_Y = -0.3,A0_on_A = 2),
#'                          register_format = TRUE)
#' x <- rtmle_init(intervals = 2,name_id = "id",name_outcome = "Y",
#'                 name_competing = "Dead",name_censoring = "Censored",
#'                 censored_label = "censored")
#' x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
#' add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
#' x <- long_to_wide(x,intervals = seq(0,2000,30.45*6))
#' protocol(x) <- list(name = "Always_A",treatment_variables = "A",intervention = 1)
#' protocol(x) <- list(name = "Never_A",treatment_variables = "A",intervention = 0)
#' prepare_data(x) <- list(reset = TRUE)
#' target(x) <- list(name = "Outcome_risk",strategy = "additive",
#'                   estimator = "tmle",protocols = c("Always_A","Never_A"))
#' x <- run_rtmle(x,time_horizon=1:2)
#' summary(x)
#'
#' @export
#' @method summary rtmle
summary.rtmle <- function(object,analysis = "Main_analysis",targets,reference = NULL,digits = 1,...){
    Reference <- risk_ratio_estimate <- risk_difference_estimate <- Protocol <- Target <- Estimate <- Upper <- Lower <- Target_parameter <- Time_horizon <- x <- NULL
    if (object$continuous_outcome){
        Target_parameter_label <- c("Mean outcome difference among survivors", "Mean outcome ratio among survivors")
    } else {
        Target_parameter_label <- c("Risk_difference", "Risk_ratio")
    }
    if (!(analysis %in% names(object$estimate))){
        stop(paste0("The object does not contain an analysis called '",
                    analysis,
                    "'. It contains:\n",
                    paste0(names(object$estimate),collapse = ", "),
                    "."))
    }
    if (length(object$estimate) == 0) {
        message("The object contains no estimates. To obtain estimates You have to call run_rtmle first. See examples.")
        return(NULL)
    }
    if (missing(targets)) {
        targets <- names(object$targets)
    }
    else {
        stopifnot(all(targets %in% names(object$targets)))
    }
    outcome_scale <- function(x){x}
    if (object$continuous_outcome == FALSE) outcome_scale <- function(x){pmin(pmax(0,100*x),100)}
    sum <- do.call(rbind,lapply(targets,function(target_name){
        target <- object$targets[[target_name]]
        protocols <- target$protocols
        risk <- do.call(rbind,lapply(protocols,function(protocol_name){
            # to avoid the internal selfdetect problem we take a copy
            e <- data.table::copy(object$estimate[[analysis]][Target == target_name & Protocol == protocol_name])
            e[,"Estimate (CI_95)":= Publish::formatCI(x = outcome_scale(e$Estimate),lower = outcome_scale(e$Lower),upper = outcome_scale(e$Upper),show.x = TRUE,digits = digits)]
            e[]
        }))
        if (analysis == "Main_analysis"){
            subset_variable <- NULL
        } else{
            subset_variable <- attr(object$estimate[[analysis]],which = "variable")
        }
        # contrasting protocols
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
                    reference_estimate <- object$estimate[[analysis]][Target == target_name & Protocol == ref & Time_horizon == tp]$Estimate
                    if (analysis == "Main_analysis"){
                        analysis_levels <- 1
                        reference_IC <- list(object$IC[[target_name]][[ref]][[paste0("time_horizon_",tp)]])
                    }else{
                        ## subset_levels <- attr(object$estimate[[analysis]],which = "levels")
                        subset_IC <- attr(object$estimate[[analysis]],which = "IC")
                        analysis_levels <- names(subset_IC)
                        if (all(analysis_levels == "")) analysis_levels <- 1:length(analysis_levels)
                        reference_IC <- lapply(analysis_levels,function(level){
                            subset_IC[[level]][[target_name]][[ref]][[paste0("time_horizon_",tp)]]                        
                        })
                        names(reference_IC) <- names(subset_IC)
                    }
                    N <- NROW(reference_IC[[1]])
                    if (analysis == "Main_analysis"&&length(object$estimate$Cheap_bootstrap)>0){
                        reference_boot <- object$estimate[["Cheap_bootstrap"]][Target == target_name & Protocol == ref & Time_horizon == tp]$Bootstrap_estimate
                    }
                    this_estimate <- object$estimate[[analysis]][Target == target_name & Protocol == protocol_name & Time_horizon == tp]$Estimate
                    if (analysis == "Main_analysis"){
                        list(this_IC <- object$IC[[target_name]][[protocol_name]][[paste0("time_horizon_",tp)]])
                    }else{
                        this_IC <- lapply(analysis_levels,function(level){
                            subset_IC[[level]][[target_name]][[protocol_name]][[paste0("time_horizon_",tp)]]                        
                        })
                        names(this_IC) <- names(subset_IC)
                    }
                    # risk difference and risk ratio
                    e <- do.call(rbind,lapply(1:length(analysis_levels),function(level){
                        risk_difference_estimate <- this_estimate[[level]] - reference_estimate[[level]]
                        risk_ratio_estimate <- exp(log(this_estimate[[level]]) - log(reference_estimate[[level]]))
                        risk_difference_se <- sd(this_IC[[level]] - reference_IC[[level]])/sqrt(N)
                        risk_difference_lower <- risk_difference_estimate - qnorm(.975)*risk_difference_se
                        risk_difference_upper <- risk_difference_estimate + qnorm(.975)*risk_difference_se
                        risk_ratio_log_se <- sd(this_IC[[level]]/this_estimate[[level]] - reference_IC[[level]]/reference_estimate[[level]])/sqrt(N)
                        risk_ratio_lower <- risk_ratio_estimate*exp(-qnorm(.975)*risk_ratio_log_se)
                        risk_ratio_upper <- risk_ratio_estimate*exp(qnorm(.975)*risk_ratio_log_se)
                        e1 <- data.table(Target = rep(target_name,2),
                                         Protocol = rep(protocol_name, 2),
                                         Target_parameter=Target_parameter_label,
                                         Time_horizon = rep(tp,2),
                                         Estimator = rep(object$estimate[[analysis]][Target == target_name & Protocol == ref]$Estimator[[1]],2),
                                         Reference = rep(ref, 2),
                                         Estimate = c(risk_difference_estimate, risk_ratio_estimate),
                                         Standard_error = c(risk_difference_se, risk_ratio_log_se),
                                         Lower = c(risk_difference_lower, risk_ratio_lower),
                                         Upper = c(risk_difference_upper, risk_ratio_upper),
                                         P_value = c(2*pnorm(-abs(risk_difference_estimate/risk_difference_se)), 2*pnorm(-abs(log(risk_ratio_estimate)/risk_ratio_log_se))))

                        if (analysis == "Main_analysis" && length(object$estimate$Cheap_bootstrap)>0){
                            this_boot <- object$estimate[["Cheap_bootstrap"]][Target == target_name & Protocol == protocol_name & Time_horizon == tp]$Bootstrap_estimate
                            boot_difference <- this_boot-reference_boot
                            boot_ratio <- this_boot/reference_boot
                            cheap_scale <- attr(object$estimate[["Cheap_bootstrap"]],"cheap_scale")
                            tq <- attr(object$estimate[["Cheap_bootstrap"]],"tq")
                            cheap_variance <- mean((risk_difference_estimate-boot_difference)^2)
                            risk_difference_boot_lower <- risk_difference_estimate - tq * cheap_scale * cheap_variance
                            risk_difference_boot_upper <- risk_difference_estimate + tq * cheap_scale * cheap_variance
                            cheap_log_variance <- mean((log(risk_ratio_estimate)-log(boot_ratio))^2)
                            risk_ratio_boot_lower <- exp(log(risk_ratio_estimate) - tq * cheap_scale * cheap_log_variance)
                            risk_ratio_boot_upper <- exp(log(risk_ratio_estimate) + tq * cheap_scale * cheap_log_variance)
                            data.table::set(e1,j = "Bootstrap_lower",value = c(risk_difference_boot_lower,risk_ratio_boot_lower))
                            data.table::set(e1,j = "Bootstrap_upper",value = c(risk_difference_boot_upper,risk_ratio_boot_upper))
                            data.table::set(e1,j = "Bootstrap_standard_error",value = c(sqrt(cheap_variance),sqrt(cheap_log_variance)))
                        }
                        if (analysis != "Main_analysis"){
                            data.table::set(e1,j = subset_variable,value = analysis_levels[[level]])
                            data.table::setcolorder(e1,subset_variable)
                        }
                        e1
                    }))
                    # risk difference
                    e[Target_parameter == Target_parameter_label[[1]],"Estimate (CI_95)":= Publish::formatCI(x = Estimate,
                                                                                                             lower = Lower,
                                                                                                             upper = Upper,
                                                                                                             show.x = TRUE,
                                                                                                             digits = digits)]
                    # risk ratio
                    e[Target_parameter == Target_parameter_label[[2]],"Estimate (CI_95)":= Publish::formatCI(x = Estimate,
                                                                                                             lower = Lower,
                                                                                                             upper = Upper,
                                                                                                             show.x = TRUE,
                                                                                                             digits = digits)]
                    e[]
                }))}))
           out <- data.table::rbindlist(list(risk,contrast),
                                         use.names = TRUE,
                                         fill = TRUE)
            out[is.na(Reference),Reference := ""][]
            if (length(subset_variable)>0){
                data.table::setkeyv(out,c(subset_variable,"Target_parameter","Protocol"))
            }else{
                data.table::setkeyv(out,cols = c("Target_parameter","Protocol"))
            }
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
