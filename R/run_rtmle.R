### run_rtmle.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  1 2024 (09:11)
## Version:
## Last-Updated: sep 11 2026 (09:27) 
##           By: Thomas Alexander Gerds
##     Update #: 643
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
#' Runs the analysis defined by \code{\link{rtmle_init}},
#' \code{\link{regime}}, \code{\link{target}}, and
#' \code{\link{model_formula}}.
#'
#' @param x Object of class \code{"rtmle"}.
#' @param targets Selection of targets to analyze. If missing, all targets in
#'     \code{x$targets} are analyzed.
#' @param learner A function, function name, or learner specification used to
#'   fit nuisance-parameter models. Must be one of:
#'   \itemize{
#'     \item A single string giving the name of a learner function
#'       (e.g., \code{"learn_glmnet"}) or the function itself.
#'     \item A list consisting of the elements:
#'       \itemize{
#'         \item \code{folds}: a character string giving the name of the
#'           learner function (e.g., \code{"learn_glmnet"}) or the function itself.
#'         \item \code{learners}: a named list of learner specifications,
#'           each of which can either be:
#'           \itemize{
#'             \item A single string.
#'             \item A list containing \code{fun} and parameters to be
#'               passed to the learner function.
#'           }
#'       }
#'   }
#' @param estimator Character specifying the estimator: either
#'     \code{"tmle"} or \code{"g-formula"}.
#' @param time_horizon The time horizon at which to calculate
#'     risks. If this is a vector, the analysis is performed for
#'     each element.
#' @param refit Logical. If \code{TRUE}, ignore any propensity score
#'     and censoring models learned in previous calls to this
#'     function. This may be useful to save computation time. The default
#'     is \code{TRUE}.
#' @param seed Seed used for cross-fitting.
#' @param subsets A list structure for subset analyses. Each element
#'     is a list requiring a label, to name the subset, and a
#'     subset of the variable \code{x$names$id} in the data
#'     \code{x$prepared_data} to identify the subset. The results of
#'     the subset analysis are stored in
#'     \code{x$estimate[[subsets[[label]]]]}. An optional element of
#'     each subset list is called \code{append}. If a result with the same
#'     label already exists, estimates are appended by default; set
#'     \code{append = FALSE} to replace the existing result. Appending may be
#'     used for stratified analyses, seed-dependence studies (Monte Carlo
#'     error), and bootstrap analyses. See
#'     examples.
#' @param keep_influence Logical. If \code{TRUE}, store the estimated
#'     influence function of the estimator in the object.  Currently
#'     this argument is used only when argument \code{subsets} is also
#'     specified.
#' @param save_fitted_objects Logical. If \code{FALSE}, store the learner
#'   summary returned as element \code{fit}. If \code{TRUE}, store the full
#'   fitted object returned as element \code{object}. Changing this setting
#'   forces nuisance-parameter models to be refitted.
#' @param progressbar Logical. If \code{TRUE}, show progress of the loops that
#'   fit the nuisance-parameter models.
#' @param verbose Logical. If \code{FALSE} suppress all
#'     messages. \code{FALSE} is the default.
#' @param ... Further arguments can change tuning parameters. For example, use
#'   \code{weight_truncation = c(0.01, 0.99)} to apply weight truncation.
#' @return The modified object contains the fitted nuisance parameter
#'     models and the estimate of the target parameter.
#' @details If the learner specification differs from the learner used for a
#'   previous run on the same object, learner-dependent estimates, fitted
#'   objects, diagnostics, and derived intervention probabilities are cleared
#'   before the new analysis is run.
#' @seealso \code{\link{rtmle_init}}, \code{\link{prepare_rtmle_data}},
#'   \code{\link{regime}}, \code{\link{target}}, \code{\link{model_formula}},
#'   \code{\link{learn_glm}}, \code{\link{learn_glmnet}},
#'   \code{\link{superlearn}}, \code{\link{summary.rtmle}}
#' @author Thomas A Gerds \email{tag@@biostat.ku.dk}
#' @examples
#' tau <- 3
#' data(simulated_cohort)
#' ld <- register_format(simulated_cohort)
#' x <- rtmle_init(time_grid = seq(0,20,4),name_id = "id",
#'                 name_outcome = "stroke",
#'                 name_competing = "death",
#'                 name_censoring = "dropout",censored_label = "censored")
#' x <- add_long_data(x,
#'                    outcome_data=ld$timevar_data$stroke[!duplicated(id)],
#'                    censored_data=ld$timevar_data$dropout,
#'                    competing_data=ld$timevar_data$death,
#'                    timevar_data=ld$timevar_data[c("bleeding","changeSBP","A","B")])
#' x <- add_baseline_data(x,data=ld$baseline_data)
#' x <- discretize_data(x,start_followup_date=0)
#' x <- prepare_rtmle_data(x)
#' x <- regime(x,name = "Always_A",
#'                     intervention = data.frame(time=x$intervention_nodes,
#'                                                    "A" = factor("1",levels = c("0","1"))))
#' x <- regime(x,name = "Never_A",
#'                     intervention = data.frame(time=x$intervention_nodes,
#'                                               "A" = factor("0",levels = c("0","1"))))
#' x <- target(x,name = "Outcome_risk",
#'                   estimator = "tmle",
#'                   regimes = c("Always_A","Never_A"))
#' x <- model_formula(x)
#' # default is undersmoothing which means: take the smallest penalty
#' # where the model still converges
#' x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau)
#' # with weight truncation
#' xw <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau,
#' weight_truncation=c(0.4,0.6))
#' # can also use lambda.min or lambda.1se
#' \dontrun{
#'     x <- run_rtmle(x,learner = list(name = "glmnet_min",
#'                                    fun="learn_glmnet",
#'                                    selector="min"),
#'                   time_horizon = tau)
#' summary(x)
#' }
#' \dontrun{
#' # Super learner combining elastic-net glmnet and ranger
#' x <- run_rtmle(
#'     x,
#'     learner = list(
#'         folds = 10,
#'         ensemble_method = "ipa",
#'         learners = list(
#'             glmnet = list(fun = "learn_glmnet",
#'                           selector = "min",
#'                           alpha = 0.5),
#'             ranger = list(fun = "learn_ranger",
#'                           num.trees = 50,
#'                           min.node.size = 10)
#'         )
#'     ),
#'     time_horizon = tau
#' )
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
                      save_fitted_objects = FALSE,
                      progressbar = 0,
                      verbose = FALSE,
                      ...){
    if (length(x$models) == 0) {
        stop(paste0("\nTODO: Use the function 'model_formula' to initialize the formula for the nuisance parameter models."))
    }
    time <- label <- level <- NULL
    if (length(x$targets) == 0) stop("Object contains no targets. You can add one with the function rtmle::target")
    dot_args <- list(...)
    x$progressbar <- progressbar
    ## x$runtime_start <- Sys.time()
    new_tuning_parms <- intersect(names(dot_args), names(x$tuning_parameters))
    if (length(new_tuning_parms)>0){
        for (ntp in new_tuning_parms) {
            x$tuning_parameters[[ntp]] <- dot_args[[ntp]]
        }
    }
    learners <- parse_learners(learner)
    learner_changed <- length(x$learner) > 0L &&
        !identical(learners, x$learner)
    if (learner_changed) {
        # Keep model formulas and regime definitions, but discard anything
        # that was fitted or calculated under the previous learner.
        clear_fitted_components <- function(value) {
            if (!is.list(value)) return(value)
            if (!is.null(names(value))) {
                value <- value[!names(value) %chin% c("fit", "fit_summary")]
            }
            lapply(value, clear_fitted_components)
        }
        x$models <- clear_fitted_components(x$models)
        x$estimate <- NULL
        x$IC <- NULL
        x$sequential_outcome_regression <- NULL
        x$run_time_horizons <- NULL
        x$diagnostics <- NULL
        x$learner <- NULL
        x$unparsed_learner <- NULL
        x$save_fitted_objects <- FALSE
        for (regime_name in names(x$regimes)) {
            x$regimes[[regime_name]]$intervention_probs <- NULL
            x$regimes[[regime_name]]$cumulative_intervention_probs <- NULL
            x$regimes[[regime_name]]$ipw_last_nodes <- NULL
            x$regimes[[regime_name]]$intervention_last_nodes <- NULL
        }
        if (length(x$prepared_data) > 0L &&
            "rtmle_predicted_outcome" %in% names(x$prepared_data)) {
            data.table::set(x$prepared_data,
                            j = "rtmle_predicted_outcome",
                            value = NULL)
        }
    }
    if (length(subsets)>0){
        refit <- TRUE
        for (sub in subsets){
            stopifnot(is.character(sub$label[[1]]))
            xs <- data.table::copy(x[c("targets","names","time_grid","time_grid_scale","time_grid_labels","regimes","models","intervention_nodes","tuning_parameters")])
            # Use a data.table as the join input.  A vector i together with
            # `on` was accepted by older data.table releases but is rejected
            # by current versions.  Keeping the index also preserves the
            # order (and possible repeated ids) of bootstrap subsets.
            subset_ids <- data.table::data.table(subset_id = sub$id)
            data.table::setnames(subset_ids, "subset_id", x$names$id)
            source_ids <- x$prepared_data[[x$names$id]]
            subset_rows <- match(sub$id, source_ids)
            subset_rows <- subset_rows[!is.na(subset_rows)]
            for (pp in names(xs$regimes)){
                xs$regimes[[pp]]$intervention_match <-
                    xs$regimes[[pp]]$intervention_match[subset_rows,
                                                          , drop = FALSE]
                xs$regimes[[pp]]$intervention_probs <- NULL
                xs$regimes[[pp]]$cumulative_intervention_probs <- NULL
            }
            # Allow for bootstrap with replacement while retaining the
            # original subject-row order.
            xs$prepared_data <- x$prepared_data[
                subset_ids,
                on = x$names$id,
                nomatch = 0L,
                allow.cartesian = TRUE
            ]
            if (NROW(xs$prepared_data) == 0) {
                stop(paste0("No data in subset: ", sub$label[[1]]))
            }
            xs$followup <- x$followup[
                subset_ids,
                on = x$names$id,
                nomatch = 0L,
                allow.cartesian = TRUE
            ]
            xs <- run_rtmle(xs,
                            targets = targets,
                            time_horizon = time_horizon,
                            seed = seed,
                            learner = learner,
                            refit = TRUE,
                            subsets = NULL,
                            save_fitted_objects = save_fitted_objects,
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
            replace_subset <- length(x$estimate[[sub$label[[1]]]]) == 0L ||
                (length(sub$append) > 0L && !isTRUE(sub$append[[1]]))
            if (replace_subset) {
                x$estimate[[sub$label[[1]]]] <- subset_result
            }else{
                # append results
                # cases stratified, Monte-Carlo error, bootstrap
                sub_IC <- attr(x$estimate[[sub$label[[1]]]],"IC",exact = TRUE)
                sub_variable <- attr(x$estimate[[sub$label[[1]]]],"variable",exact = TRUE)
                ## sub_level <- attr(x$estimate[[sub$label[[1]]]],"level",exact = TRUE)
                x$estimate[[sub$label[[1]]]] <- rbind(x$estimate[[sub$label[[1]]]],
                                                      subset_result,
                                                      fill = TRUE)
                data.table::setattr(x$estimate[[sub$label[[1]]]],
                                    "IC",
                                    c(sub_IC, attr(subset_result,"IC",exact = TRUE)))
                data.table::setattr(x$estimate[[sub$label[[1]]]],
                                    "variable",sub_variable)# assumed equal to attr(subset,"variable",exact = TRUE)
                ## data.table::setattr(x$estimate[[sub$label[[1]]]],
                ## "level",
                ## c(sub_level, attr(subset_result,"level",exact = TRUE)))
            }
        }
        x$unparsed_learner <- learner
        x$learner <- learners
        x$save_fitted_objects <- save_fitted_objects
        return(x)
    }else{
        #
        # check data
        #
        # Skipping the nuisance parameter models is only possible when the same
        # learner was used previously
        if (length(x$learner) == 0L || learner_changed){
            refit <- TRUE
        }
        if (!identical(isTRUE(x$save_fitted_objects),isTRUE(save_fitted_objects))){
            refit <- TRUE
        }
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
            time_horizon <- max(x$time_grid)
        } else {
            stopifnot(all(time_horizon <= max(x$time_grid) & time_horizon>0))
        }
        if (!(x$names$id%in%names(x$prepared_data)))
            stop(paste0("Cannot see id variable ",x$names$id," in x$prepared_data."))
        ## make sure that the treatment variables are factors with levels equal to
        # those specified by the regimes
        for (v in names(x$names$treatment_options)){
            v_treatment_variables <- intersect(paste0(v,"_",x$time_grid),names(x$prepared_data))
            for (v_j in v_treatment_variables){
                if (inherits(x$prepared_data[[v_j]],"factor")){
                    if (!(all.equal(levels(x$prepared_data[[v_j]]),as.character(x$names$treatment_options[[v]])))){
                        stop(paste0("The regimes specify the following treatment options (factor levels) for variable ",v,
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
            stats::setNames(lapply(x$targets[[target_name]]$regimes,function(regime_name){
                stats::setNames(lapply(1:length(time_horizon),function(th){
                    numeric(NROW(x$prepared_data))
                }),label_time_horizon)
            }),x$targets[[target_name]]$regimes)}),run_these_targets)
        # initialize estimate table
        empty_estimate <- data.table::rbindlist(lapply(run_these_targets,function(target_name){
            data.table::rbindlist(lapply(x$targets[[target_name]]$regimes,function(regime_name){
                expand.grid(Target = target_name,
                            Regime = regime_name,
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
            # initialize new targets, new regimes and new time_horizons
            e <- rbind(x$estimate[["Main_analysis"]],empty_estimate)
            e <- e[e[,.I[1],by = c("Target","Regime","Time_horizon","Estimator")]$V1]
            x$estimate[["Main_analysis"]] <- e
        }
        # for loop across regimes
        run_these_regimes <- unique(unlist(sapply(run_these_targets,function(target_name){
            x$targets[[target_name]]$regimes})))
        missing_regimes <- setdiff(run_these_regimes, names(x$regimes))
        if (length(missing_regimes) > 0L) {
            missing_by_target <- vapply(
                run_these_targets,
                function(target_name) {
                    target_missing <- intersect(
                        x$targets[[target_name]]$regimes,
                        missing_regimes
                    )
                    if (length(target_missing) == 0L) {
                        return("")
                    }
                    paste0(
                        "target '", target_name, "': ",
                        paste(target_missing, collapse = ", ")
                    )
                },
                character(1L)
            )
            missing_by_target <- missing_by_target[nzchar(missing_by_target)]
            available_regimes <- names(x$regimes)
            if (length(available_regimes) == 0L) {
                available_regimes <- "(none)"
            }
            stop(
                "Cannot run rtmle because the following regime(s) are not ",
                "defined: ", paste(missing_regimes, collapse = ", "),
                ". Referenced by ", paste(missing_by_target, collapse = "; "),
                ". Available regimes: ",
                paste(available_regimes, collapse = ", "), "."
            )
        }
        # use the first regime to store the censoring models
        if (length(x$censoring_use_regime) == 0
            || !(x$censoring_use_regime %chin% names(x$regimes))){
            x$censoring_use_regime <- run_these_regimes[[1]]
        }
        if (estimator == "tmle"){
            for (regime_name in run_these_regimes){
                #
                # G-part: fit nuisance parameter models for propensity and censoring
                #
                # when regimes are defined before data are prepared then
                # intervention_match needs to run here
                x <- intervention_match(x,regime_name = regime_name)
                x <- intervention_probabilities(x,
                                                regime_name = regime_name,
                                                max_intervention_node = max(time_horizon)-1,
                                                refit = refit,
                                                learner = learners,
                                                seed = seed,
                                                progressbar = progressbar,
                                                save_fitted_objects = save_fitted_objects)
            }
        }
        if (save_fitted_objects){
            # FIXME: go over saved censoring models and replace the objects with their summary
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
            for (regime_name in x$targets[[target_name]]$regimes){
                if (verbose[[1]]){
                    message("Current regime: ",regime_name," ... Set argument verbose = FALSE to suppress this message.")
                }
                #
                # Q-part: loop backwards in time through iterative condtional expectations
                #
                # loop across time-horizons
                for (th in time_horizon){
                    if (verbose){
                        message("Running sequential regression backwards from time ",th)
                    }
                    x <- sequential_regression(x = x,
                                               target_name = target_name,
                                               regime_name = regime_name,
                                               time_horizon = th,
                                               learner = learners,
                                               estimator = estimator,
                                               seed = seed,
                                               progressbar = progressbar,
                                               save_fitted_objects = save_fitted_objects)
                }
            }
        }
        ## store the time_horizon
        x$run_time_horizons <- unique(c(x$run_time_horizons,time_horizon))
        ## Keep the learner function used for cheap bootstrap
        x$unparsed_learner <- learner
        x$learner <- learners
        x$save_fitted_objects <- save_fitted_objects
        return(x)
    }
}

######################################################################
### run_rtmle.R ends here
