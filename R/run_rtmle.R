### run_rtmle.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  1 2024 (09:11)
## Version:
## Last-Updated: Mar 28 2025 (13:53) 
##           By: Thomas Alexander Gerds
##     Update #: 393
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
#' @param seed Seed used for cross-fitting
#' @param ... Additional arguments passed to the learner function.
#' @return The modified object contains the fitted nuisance parameter models and the estimate of the target parameter.
#' @author  Thomas A Gerds \email{tag@@biostat.ku.dk} 
#' @examples
#' # ------------------------------------------------------------------------------------------
#' # Intervening on a single treatment variable
#' # ------------------------------------------------------------------------------------------
#'
#' set.seed(17)
#' tau <- 3
#' ld <- simulate_long_data(n = 91,number_visits = 20,
#'                          beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),
#'                          register_format = TRUE)
#' x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
#' add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
#' x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
#' protocol(x) <- list(name = "Always_A",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1"))))
#' protocol(x) <- list(name = "Never_A",
#'                     intervention = data.frame("A" = factor("0",levels = c("0","1"))))
#' prepare_data(x) <- list()
#' target(x) <- list(name = "Outcome_risk",
#'                   strategy = "additive",
#'                   estimator = "tmle",
#'                   protocols = c("Always_A","Never_A"))
#' x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau)
#' summary(x)
#'
#'
#'
#' @export
run_rtmle <- function(x,
                      targets,
                      learner = "learn_glm",
                      time_horizon,
                      refit = FALSE,
                      seed = NULL,
                      ...){
    time <- value <- NULL
    if (length(x$continuous_outcome) == 0 || x$continuous_outcome == FALSE){
        x$continuous_outcome <- FALSE
    }else{
        x$continuous_outcome <- TRUE
    }
    # check data
    ## sapply(x$prepared_data,function(x)sum(is.na(x)))
    if (!(x$names$id%in%names(x$prepared_data)))
        stop(paste0("Cannot see id variable ",x$names$id," in x$prepared_data."))
    ## make sure that the treatment variables are factors with levels equal to those specified by the protocols
    ## FIXME: this could perhaps be done in prepare_data()
    for (v in names(x$names$treatment_options)){
        for (v_j in paste0(v,"_",x$times)){
            if (inherits(x$prepared_data[[v_j]],"factor")){
                if (!(all.equal(levels(x$prepared_data[[v_j]]),as.character(x$names$treatment_options[[v]])))){
                    stop(paste0("The protocols specify the following treatment options (factor levels) for variable ",v,
                                paste0(x$names$treatment_options[[v]],collapse = ","),"\nBut, the data have: ",
                                paste0(levels(x$prepared_data[[v_j]]),collapse = ",")))
                }
            }else{
                ## stop(paste0("The treatment variable ",v_j," is not a factor"))
                data.table::set(x$prepared_data,
                                j = v_j,
                                value = factor(x$prepared_data[[v_j]],levels = x$names$treatment_options[[v]]))
            }
        }
    }
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
        time_horizon <- max(x$times)
    } else {
        stopifnot(all(time_horizon <= max(x$times) & time_horizon>0))
    }
    if (length(learner)>1){
        learners <- parse_learners(learner)
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
                                            seed = seed,
                                            ...)
            #
            # Q-part: loop backwards in time through iterative condtional expectations
            #
            # initialize estimate
            if (!x$continuous_outcome){
                Target_parameter <- "Risk"
            } else {
                Target_parameter <- "Weighted mean among survivors"
            }
            x$estimate[[target_name]][[protocol_name]] <- data.table(Target = target_name,
                                                                     Protocol = protocol_name,
                                                                     Target_parameter = Target_parameter,
                                                                     Time_horizon = time_horizon,
                                                                     Estimator = x$targets[[target_name]]$estimator,
                                                                     Estimate = numeric(length(time_horizon)),
                                                                     P_value = NA)
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
                                           learner = learner,
                                           seed = seed,
                                           ...)
            }
        }
    }
    # FIXME: where did x loose its class?
    class(x) <- "rtmle"
    x
}
######################################################################
### run_rtmle.R ends here
