### run_rtmle.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  1 2024 (09:11)
## Version:
## Last-Updated: dec  2 2025 (09:52) 
##           By: Thomas Alexander Gerds
##     Update #: 562
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
#' @param targets Selection of targets to be analysed. If missing all
#'     targets in x$targets are analysed.
#' @param learner A function which is called to fit (learn) the
#'     nuisance parameter models.
#' @param estimator  Character specifying the estimator: either \code{'tmle'} or \code{'g-formula'}.
#' @param time_horizon The time horizon at which to calculate
#'     risks. If it is a vector the analysis will be performed for
#'     each element of the vector.
#' @param refit Logical. If \code{TRUE} ignore any propensity score
#'     and censoring models learned in previous calls to this
#'     function. This may be useful to save computation time. Default is \code{TRUE}.
#' @param seed Seed used for cross-fitting
#' @param subsets A list structure for subset analyses. Each element is a list 
#'                which requires a label, to name the subset, and a subset of the 
#'                variable \code{x$names$id} in the data \code{x$prepared_data} to
#'                identify the subset. The results of the subset analysis are stored
#'                in \code{x$estimate[[subsets[[label]]}. An optional element of each
#'                subset-list is called \code{append}
#'                which should be logical: if \code{TRUE} append the estimates to the existing
#'                estimates with rbind. This may be used for stratified analyses, to study seed dependence (Monte-Carlo error)
#'                and bootstrap. See examples.
#' @param keep_influence Logical: if \code{TRUE} store the estimated influence function of the estimator in the object.
#'                       Currently this argument is only used when argument \code{subsets} is also specified.
#' @param verbose Logical. If \code{FALSE} suppress all messages. \code{FALSE} is the default.
#' @return The modified object contains the fitted nuisance parameter
#'     models and the estimate of the target parameter.
#' @author Thomas A Gerds \email{tag@@biostat.ku.dk}
#' @examples
#' # ------------------------------------------------------------------------------------------
#' # Intervening on a single treatment variable
#' # ------------------------------------------------------------------------------------------
#'
#' set.seed(17)
#' tau <- 3
#' ld <- simulate_long_data(n = 391,number_visits = 20,
#'                          beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),
#'                          register_format = TRUE)
#' x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- add_long_data(x,
#'                    outcome_data=ld$outcome_data,
#'                    censored_data=ld$censored_data,
#'                    competing_data=ld$competing_data,
#'                    timevar_data=ld$timevar_data)
#' x <- add_baseline_data(x,data=ld$baseline_data)
#' x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
#' x <- protocol(x,name = "Always_A",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1"))))
#' x <- protocol(x,name = "Never_A",
#'                     intervention = data.frame("A" = factor("0",levels = c("0","1"))))
#' x <- prepare_data(x)
#' x <- target(x,name = "Outcome_risk",
#'                   estimator = "tmle",
#'                   protocols = c("Always_A","Never_A"))
#' x <- model_formula(x)
#' # default is undersmoothing which means: take the smallest penalty
#' # where the model still converges
#' x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau)
#' # can also use lambda.min or lambda.1se
#' \dontrun{
#' x <- run_rtmle(x,learner = list("glmnet_cv"=list(learner_fun="learn_glmnet",
#'                                 selector="min")),
#'               time_horizon = tau)
#' summary(x)
#' }
#' \dontrun{
#' # stratified analyses
#' x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = tau,
#'                verbose=FALSE,
#'                subsets=list(list(label="Sex",variable="Sex",
#'                                  level="Female",id=x$prepared_data[sex==0,id]),
#'                        list(label="Sex",variable="Sex",
#'                                  level="Male",id=x$prepared_data[sex==1,id])))
#' }
#'
#'
#'
#' @export
run_rtmle <- function(x,
                      targets,
                      learner = "learn_glm",
                      estimator = "tmle",
                      time_horizon,
                      refit = TRUE,
                      seed = NULL,
                      subsets = NULL,
                      keep_influence = TRUE,
                      verbose = FALSE){
    time <- label <- level <- NULL
    if (length(x$targets) == 0) stop("Object contains no targets. You can add one with the function rtmle::target")
    ## FIXME: need to stop or adapt if learner="glmnet" instead of "learn_glmnet"
    if (length(subsets)>0){
        for (sub in subsets){
            stopifnot(is.character(sub$label[[1]]))
            ## FIXME: is this check useful? stopifnot(all(sub$id %in% x$prepared_data[[x$names$id]]))
            xs <- data.table::copy(x[c("targets","names","times","protocols","models")])
            xs$prepared_data <- x$prepared_data[x$prepared_data[[x$names$id]] %in% sub$id]
            if (NROW(xs$prepared_data) == 0) stop(paste0("No data in subset: ",label))
            # FIXME: should check if the number at risk at the maximal time_horizon is not zero
            #        and that there is still variation in the outcome at that time
            xs$followup <- x$followup[x$followup[[x$names$id]] %in% sub$id]
            xs <- run_rtmle(xs,
                            targets = targets,
                            time_horizon = time_horizon,
                            seed = seed,
                            learner = learner,
                            refit = TRUE,
                            subsets = NULL,
                            verbose = verbose)
            subset_result <- xs$estimate[["Main_analysis"]]
            # add the subset identifying information, such as age="40-60"
            if (length(sub$variable)>0){
                if (length(sub$level)>0) v <- sub$level else v = ""
                addthis <- data.table::data.table(v)
                data.table::setnames(addthis,sub$variable[[1]])
                subset_result <- cbind(addthis,subset_result)
                data.table::setattr(subset_result,"variable",sub$variable)
                ## data.table::setattr(subset_result,"level",sub$level)
            }
            # add the fitted nuisance parameter models
            for (m in names(x$models)){
                x$models[[m]][[sub$label]] <- xs$models[[m]][["fit"]]
            }
            # add the estimate of the influence function
            if (keep_influence) {
                vic <- list(xs$IC)
                names(vic) <- if (length(sub$level)>0) sub$level else ""
                data.table::setattr(subset_result,"IC",vic)
            }
            # set or replace existing results
            if (length(x$estimate[[sub$label[[1]]]]) == 0 ||
                (length(sub$append) == 0) ||
                sub$append[[1]] == FALSE){
                x$estimate[[sub$label[[1]]]] <- subset_result
            }else{
                # append results
                # cases stratified, Monte-Carlo error, bootstrap
                sub_IC <- attr(x$estimate[[sub$label[[1]]]],"IC")
                sub_variable <- attr(x$estimate[[sub$label[[1]]]],"variable")
                ## sub_level <- attr(x$estimate[[sub$label[[1]]]],"level")
                x$estimate[[sub$label[[1]]]] <- rbind(x$estimate[[sub$label[[1]]]],
                                                      subset_result,
                                                      fill = TRUE)
                data.table::setattr(x$estimate[[sub$label[[1]]]],
                                    "IC",
                                    c(sub_IC, attr(subset_result,"IC")))
                data.table::setattr(x$estimate[[sub$label[[1]]]],
                                    "variable",sub_variable)# assumed equal to attr(subset,"variable")
                ## data.table::setattr(x$estimate[[sub$label[[1]]]],
                ## "level",
                ## c(sub_level, attr(subset_result,"level")))
            }
        }
        return(x)
    }else{
        # check data
        Target_parameter <- "Risk"
        available_targets <- names(x$targets)
        if (!missing(targets) && length(targets)>0) {
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
        learners <- parse_learners(learner)
        ## sapply(x$prepared_data,function(x)sum(is.na(x)))
        if (!(x$names$id%in%names(x$prepared_data)))
            stop(paste0("Cannot see id variable ",x$names$id," in x$prepared_data."))
        ## make sure that the treatment variables are factors with levels equal to those specified by the protocols
        ## FIXME: this could perhaps be done in prepare_data()
        for (v in names(x$names$treatment_options)){
            v_treatment_variables <- intersect(paste0(v,"_",x$times),names(x$prepared_data))
            for (v_j in v_treatment_variables){
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
                                    value = factor(x$prepared_data[[v_j]],
                                                   levels = x$names$treatment_options[[v]]))
                }
            }
        }
        # initialize object to receive estimates and IC 
        label_time_horizon <- paste0("time_horizon_",time_horizon)
        x$sequential_outcome_regression <- vector(mode = "list",length(run_these_targets))
        names(x$sequential_outcome_regression) = run_these_targets
        # initialize influence curve vector
        x$IC <- stats::setNames(lapply(run_these_targets,function(target_name){
            stats::setNames(lapply(x$targets[[target_name]]$protocols,function(protocol_name){
                stats::setNames(lapply(1:length(time_horizon),function(th){
                    numeric(NROW(x$prepared_data))
                }),label_time_horizon)
            }),x$targets[[target_name]]$protocols)}),run_these_targets)
        # initialize estimate table
        empty_estimate <- data.table::rbindlist(lapply(run_these_targets,function(target_name){
            data.table::rbindlist(lapply(x$targets[[target_name]]$protocols,function(protocol_name){
                expand.grid(Target = target_name,
                            Protocol = protocol_name,
                            Target_parameter = Target_parameter,
                            Time_horizon = time_horizon,
                            Estimator = estimator,
                            Estimate = numeric(1),
                            P_value = 1,
                            Standard_error = numeric(1),
                            Lower = numeric(1),
                            Upper = numeric(1))
            }))}))
        if (length(x$estimate[["Main_analysis"]]) == 0){
            x$estimate[["Main_analysis"]] <- empty_estimate
        }else{
            # initialize new targets, new protocols and new time_horizons
            e <- rbind(x$estimate[["Main_analysis"]],empty_estimate)
            e <- e[e[,.I[1],by = c("Target","Protocol","Time_horizon")]$V1]
            x$estimate[["Main_analysis"]] <- e
        }
        # for loop across targets
        for (target_name in run_these_targets){
            if (verbose[[1]]){
                message("Running target: ",
                        target_name,
                        "... Set argument verbose = FALSE to suppress this message.")
            }
            x$sequential_outcome_regression[[target_name]] <- vector(mode = "list",3)
            names(x$sequential_outcome_regression[[target_name]]) <- c("predicted_values","fit","intervened_data")
            for (protocol_name in x$targets[[target_name]]$protocols){
                if (verbose[[1]]){
                    message("Current protocol: ",protocol_name," ... Set argument verbose = FALSE to suppress this message.")
                }
                #
                # G-part: fit nuisance parameter models for propensity and censoring
                #
                x <- intervention_probabilities(x,
                                                protocol_name = protocol_name,
                                                refit = refit,
                                                learner = learners,
                                                time_horizon = time_horizon,
                                                seed = seed)
                #
                # Q-part: loop backwards in time through iterative condtional expectations
                #
                # loop across time-horizons
                for (th in time_horizon){
                    x <- sequential_regression(x = x,
                                               target_name = target_name,
                                               protocol_name = protocol_name,
                                               time_horizon = th,
                                               learner = learners,
                                               estimator = estimator,
                                               seed = seed)
                    # FIXME: why does x loose its class in this loop?
                    class(x) <- "rtmle"
                }
            }
        }
        ## Keep the learner function used for cheap bootstrap
        x$learner <- learner
        return(x)
    }
}

######################################################################
### run_rtmle.R ends here
