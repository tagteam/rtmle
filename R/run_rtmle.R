### run_rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  1 2024 (09:11) 
## Version: 
## Last-Updated: Nov  3 2024 (14:33) 
##           By: Thomas Alexander Gerds
##     Update #: 359
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Sequential regression with TMLE update step for discretized follow-up data
#'
#' This function runs the analysis defined in previous steps.
#' @param x object of class \code{rtmle} 
#' @param targets Selection of targets to be analysed. If missing all targets in x$targets are analysed.
#' @param learner A function which is called to fit (learn) the nuisance parameter models.
#' @param time_horizon The time horizon at which to calculate risks. If it is a vector the analysis will be performed for each element of the vector. 
#' @param refit Logical. If \code{TRUE} ignore any propensity score and censoring models learned in previous calls to this function. Default is \code{FALSE}.
#' @param ... Additional arguments passed to the learner function.
#' @export
run_rtmle <- function(x,
                      targets,
                      learner = "learn_glm",
                      time_horizon,
                      refit = FALSE,
                      ...){
    time <- value <- NULL
    # check data 
    ## sapply(x$prepared_data,function(x)sum(is.na(x)))
    requireNamespace("riskRegression")
    # for loop across targets
    available_targets <- names(x$targets)
    if (!missing(targets)) {
        if (!(all(targets %in% available_targets)))
            stop(paste0("Requested targets: \n",paste(targets,collapse = ","),"\n\navailable targets:\n",paste(available_targets,collapse = ",")))
        run_these_targets <- intersect(targets,available_targets)
    }else{
        run_these_targets <- available_targets
    }
    if (missing(time_horizon)) {
        time_horizon <- max(x$time)
    } else {
        stopifnot(all(time_horizon <= max(x$time) & time_horizon>0))
    }
    if (length(learner)>1){
        for (this_learner in 1:length(learner)){
            ## if no learner function is specified then assume that the name
            ## is a function
            if (is.character(name_is_fun <- learner[[this_learner]])){
                learner[[this_learner]] <- list(learner_fun = name_is_fun)
                names(learner)[[this_learner]] <- name_is_fun
            }else{
                if (length(learner[[this_learner]][["learner_fun"]]) == 0)
                    learner[[this_learner]]$learner_fun <- names(learner)[[this_learner]]
            }
        }
        ## make unique names
        names(learner) <- make.names(names(learner),unique = TRUE)
    }
    for (target_name in run_these_targets){
        message("Running target: ",target_name)
        x$sequential_outcome_regression[[target_name]] <- vector(mode = "list",3)
        names(x$sequential_outcome_regression[[target_name]]) <- c("predicted_values","fit","intervened_data")
        ## protocols <- names(x$protocols)
        protocols <- x$targets[[target_name]]$protocols
        for (protocol_name in protocols){
            message("Current protocol: ",protocol_name)
            #
            # G-part: fit nuisance parameter models for propensity and censoring
            #
            x <- intervention_probabilities(x,
                                            protocol_name = protocol_name,
                                            refit = refit,
                                            learner = learner,
                                            time_horizon = time_horizon,
                                            ...)
            # 
            # Q-part: loop backwards in time through iterative condtional expectations
            #
            # initialize estimate
            x$estimate[[target_name]][[protocol_name]] <- data.table(Target = target_name,
                                                                     Protocol = protocol_name,
                                                                     "Time_horizon" = time_horizon,
                                                                     Estimator = x$targets[[target_name]]$estimator,
                                                                     Estimate = numeric(length(time_horizon)))
            label_time_horizon <- paste0("time_horizon_",time_horizon)
            x$sequential_outcome_regression <- vector(mode = "list",length(run_these_targets))
            names(x$sequential_outcome_regression) = run_these_targets
            # initialize influence curve vector
            for (th in 1:length(time_horizon)){
                # FIXME: sample size
                x$IC[[target_name]][[protocol_name]][[label_time_horizon[[th]]]] <- numeric(NROW(x$prepared_data))
            }
            # loop across time-horizons
            for (th in time_horizon){
                x <- sequential_regression(x = x,
                                           target_name = target_name,
                                           protocol_name = protocol_name,
                                           time_horizon = th,
                                           learner = learner,...)
            }
        }
    }
    # FIXME: where did x loose its class?
    class(x) <- "rtmle"
    x
}
######################################################################
### run_rtmle.R ends here
