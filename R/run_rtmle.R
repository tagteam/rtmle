### run_rtmle.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  1 2024 (09:11)
## Version:
## Last-Updated: Nov 27 2026 (16:05) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 668
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
#' Sequential regression for discrete-time or continuous-time follow-up data
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
#' @param estimator Character specifying the estimator. Discrete-time analyses
#'     support \code{"tmle"} and \code{"g-formula"}. Continuous-time analyses
#'     support \code{"tmle"}, \code{"one-step"} for debiased ICE-IPCW, and
#'     \code{"sequential_regression"} for the ICE-IPCW plug-in estimate.
#' @param time_horizon The time horizon at which to calculate
#'     risks. If this is a vector, the analysis is performed for
#'     each element. These are indices into the registered time grid, starting
#'     at 1; by default the last node is used.
#' @param refit Logical. If \code{TRUE}, ignore any propensity score
#'     and censoring models learned in previous calls to this
#'     function. This may be useful to save computation time. The default
#'     is \code{TRUE}.
#' @param seed Seed used for cross-fitting.
#' @param save_fitted_objects Logical. If \code{FALSE}, store the learner
#'   summary returned as element \code{fit}. If \code{TRUE}, store the full
#'   fitted object returned as element \code{object}. Changing this setting
#'   forces nuisance-parameter models to be refitted.
#' @param progressbar Logical. If \code{TRUE}, show progress of the loops that
#'   fit the nuisance-parameter models.
#' @param verbose Logical. If \code{TRUE}, show progress messages.
#'     \code{FALSE} is the default.
#' @param ... Further arguments can change tuning parameters. For example, use
#'   \code{weight_truncation = c(0.01, 0.99)} to apply weight truncation.
#' @return The modified object contains the fitted nuisance parameter
#'     models and the estimate of the target parameter.
#' @details If the learner specification differs from the learner used for a
#'   previous run on the same object, learner-dependent estimates, fitted
#'   objects, diagnostics, and derived intervention probabilities are cleared
#'   before the new analysis is run.
#'
#' For continuous-time analyses, initialize with
#' \code{rtmle_init(continuous_time = TRUE)} and supply event histories with
#' \code{add_continuous_data}, using the separate tables accepted by
#' \code{add_long_data}. Add baseline predictors with \code{add_baseline_data}.
#' Continue with \code{prepare_rtmle_data},
#' \code{regime}, \code{target}, \code{model_formula}, and this function.
#' Entries in \code{time_grid} are expressed in the original time units;
#' \code{time_horizon} selects their node indices as in discrete time.
#' Static and dynamic interventions on binary treatment \code{A} are
#' supported at baseline and at each treatment opportunity. The default is
#' \code{conservative = TRUE}. With \code{estimator = "one-step"},
#' \code{conservative = FALSE} adds an experimental censoring-martingale
#' correction to the estimate and to each participant's influence function.
#' The default outcome learner is \code{learn_oipcw}, which fits IPCW
#' pseudo-outcomes that may exceed one. Use \code{learner = "ipcw_glm"} for
#' a quasi-binomial GLM with the unweighted pseudo-outcome as the response
#' and IPCW as observation weights. The default link is logit. For another
#' link, use a specification such as
#' \code{list(name = "ipcw_probit", fun = "ipcw_glm", link = "probit")}.
#' Specify \code{ipcw_glm} as a single outcome or grid learner.
#'
#' Continuous-time options passed through \code{...} are:
#' \itemize{
#'   \item \code{propensity_learner}: a learner specification for treatment in
#'     \code{"tmle"} and \code{"one-step"} analyses. It defaults to the supplied
#'     \code{learner}, or \code{"learn_glm"} if \code{learner} is omitted or set
#'     to \code{ipcw_glm}.
#'   \item \code{censoring_learner}: a single survival learner specification,
#'     defaulting to \code{"learn_coxph"}. Include \code{marginal_censoring}
#'     in the learner list to choose a model depending only on baseline covariates (\code{TRUE}, the
#'     default) or event-history models (\code{FALSE}). For a marginal model,
#'     \code{marginal_censoring_formula} provide a
#'     one-sided formula; Cox models may use \code{strata()}. Both model types
#'     support \code{"learn_coxph"}, \code{"learn_weibull"}, and custom learners.
#'     Marginal custom learners must accept \code{reuse_fit}, \code{time_variable},
#'     and \code{time_ref}, returning interval survival in \code{predicted_values}.
#'     Alternatively, \code{time_varying = TRUE} allows for a time-varying Cox
#'     model. Use \code{learn_coxph}, \code{learn_aalen}, or
#'     \code{learn_cox_aalen}, and optionally a one-sided
#'     \code{censoring_formula}, for example \code{~ A + L + age} for Cox/Aalen
#'     or \code{~ L + prop(A) + prop(age)} for Cox-Aalen. The default uses current
#'     treatment/covariates and fixed baseline predictors. Initial \code{*_0}
#'     values can be requested explicitly. This mode cannot be combined with
#'     \code{marginal_censoring = TRUE}. 
#'     The workflow constructs its interval data from \code{add_continuous_data}.
#'     The Aalen learners require the optional \code{timereg} package.
#'     See the continuous-time vignette for examples.
#'   \item \code{event_learner}: survival learner for the time to the next
#'     full-data event in the non-conservative analysis. The default is
#'     \code{"learn_coxph"}. Each event index uses elapsed time and the
#'     history through the preceding event. Treatment and covariate updates,
#'     outcomes and competing events count as failures. Dropout and the
#'     administrative end of follow-up right-censor the event time.
#'     The baseline-only option is unavailable for this model. Custom survival
#'     learners used for the correction must predict from a saved fit at
#'     arbitrary endpoints and honor \code{left_limit} at hazard jumps.
#'   \item \code{grid_learner}: outcome learner for interior grid times.
#'     Defaults to the outcome learner used at the analysis horizon. With
#'     \code{ipcw_glm}, each grid regression uses weights for its own cutoff.
#'   \item \code{integration_grid_size}: number of equal subintervals from
#'     zero to each requested horizon for the non-conservative correction
#'     (default 100, minimum 2). Uses additional pseudo-outcome at
#'     prior time horizons.
#'   \item \code{lag}: number of previous event indices
#'     to consider. The baseline covariates are always
#'     included.
#'   \item \code{gbound}: an optional lower bound in \eqn{(0,1]} for each
#'     treatment probability under the regime.
#'   \item \code{stratify_sequential_regression}: stratify on adherent histories
#'     (default \code{FALSE}).
#'   \item \code{intervene_all_sequential_regression}: in the sequential regression steps
#'     whether to predict using only by changing the current treatment values
#'     or by changing all.
#'     This option is currently unavailable for dynamic continuous-time regimes.
#'   \item \code{intervention_hook}: an optional function \code{function(data, k)}
#'     updating treatment-dependent predictors after intervention. It must
#'     preserve subject identifiers and row order.
#'   \item \code{exclude_latest_covariate}: covariate names whose current values
#'     should be omitted from treatment models.
#' }
#' Changing continuous-time analysis settings clears previous estimates so
#' stored influence curves and estimates refer to the same analysis.
#' Subsetting/stratifying continuous-time analyses
#' and ensemble learners for censoring are not supported. The continuous-time
#' \code{"sequential_regression"} option uses outcome regressions with inverse
#' probability of censoring weights when censoring is present.
#' It reports the ICE-IPCW plug-in estimate.
#' \code{summary} and risk plots use the usual rtmle methods. Outcome model
#' summaries are stored in \code{x$sequential_outcome_regression}; treatment and
#' censoring summaries are stored in \code{x$continuous_nuisance}.
#' @seealso \code{\link{rtmle_init}}, \code{\link{prepare_rtmle_data}},
#'   \code{\link{regime}}, \code{\link{target}}, \code{\link{model_formula}},
#'   \code{\link{learn_glm}}, \code{\link{learn_glmnet}},
#'   \code{\link{superlearn}}, \code{\link{summary.rtmle}},
#'   \code{\link{subset.rtmle}}
#' @seealso \code{\link{add_continuous_data}}, \code{\link{learn_oipcw}}, \code{\link{learn_ipcw_glm}},
#'   \code{\link{learn_coxph}}
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
#' @export
run_rtmle <- function(x,
                      targets,
                      learner = "learn_glm",
                      estimator = "tmle",
                      time_horizon,
                      refit = TRUE,
                      seed = NULL,
                      save_fitted_objects = FALSE,
                      progressbar = 0,
                      verbose = FALSE,
                      ...){
    if (isTRUE(x$continuous_time)) {
        return(run_rtmle_continuous(x,
                                    targets = if (missing(targets)) NULL else targets,
                                    learner = if (missing(learner)) NULL else learner,
                                    estimator = estimator,
                                    time_horizon = if (missing(time_horizon)) NULL else time_horizon,
                                    refit = refit, seed = seed,
                                    save_fitted_objects = save_fitted_objects,
                                    progressbar = progressbar, verbose = verbose, ...))
    }
    if (length(x$models) == 0) {
        stop(paste0("\nTODO: Use the function 'model_formula' to initialize the formula for the nuisance parameter models."))
    }
    if (length(x$targets) == 0) stop("Object contains no targets. You can add one with the function rtmle::target")
    dot_args <- list(...)
    if (any(c("subsets", "keep_influence") %in% names(dot_args))) {
        stop("Use subset() to select participants before calling run_rtmle().")
    }
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
        x <- clear_rtmle_fit(x)
    }
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

######################################################################
### run_rtmle.R ends here
