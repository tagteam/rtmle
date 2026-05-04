### superlearn.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 31 2024 (07:29) 
## Version: 
## Last-Updated: apr 29 2026 (07:25) 
##           By: Thomas Alexander Gerds
##     Update #: 297
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Super learn nuisance-parameter models
##'
##' Uses cross-fitting to produce level-one data for learners used by
##' \code{\link{run_rtmle}}. The predicted probability from each learner is
##' stored as one column and compared with the outcome. The discrete super
##' learner chooses the learner with the lowest quadratic prediction error; the
##' ensemble super learner combines the candidate learners.
##' 
##' @title Super learning for rtmle
##' @param folds Number of folds.
##' @param seed Random seed.
##' @param learners List of learners. Each learner can either be a
##'     character string or a list. If the learner is a character
##'     string it must be the name of a learner function, such as
##'     \code{\link{learn_glm}}, \code{\link{learn_glmnet}},
##'     \code{\link{learn_ranger}}, or \code{\link{learn_xgboost}}. If the
##'     learner is a list, either the name of the learner is a learner function
##'     or the list has an element \code{fun} naming a learner function. The
##'     other elements of the list are passed
##'     as additional arguments to the learner function.
##' @param parse_learners Logical. If \code{TRUE}, parse learner
##'     specifications with \code{\link{parse_learners}}.
##' @param character_formula A formula (passed as a character string!)
##'     to parse the outcome and the predictor variables.
##' @param outcome_variable The name of the outcome variable (same as
##'     \code{all.vars(formula(character_formula))[[1]]}), provided to
##'     avoid overhead in the parsing of the formula.
##' @param outcome_variable_name For internal use when called from
##'   \code{\link{run_rtmle}} via \code{fitter()}.
##' @param id_variable Name of the subject identifier variable.
##' @param data Data used for learning.
##' @param intervened_data Data in which all intervention variables have already
##'     been set according to the intervention protocol.
##' @param ensemble_method How to combine learners. Implemented methods are
##'   non-linear least squares (\code{"nnls"}), index-of-prediction-accuracy
##'   weighting (\code{"ipa"}), and \code{"discrete"}, which selects the learner
##'   with the lowest Brier score. If all models have zero weight, the average
##'   predicted value from the learning set is used.
##' @param diagnostics For internal use when called from \code{\link{run_rtmle}}
##'   via \code{fitter()}, which is called from
##'   \code{intervention_probabilities()} and \code{sequential_regression()}.
##' @param ... Not used.
##' @return A list whose first element, \code{predicted_values}, is the
##'   prediction vector from the super learner. Element \code{fit} contains the
##'   ensemble weights, and element \code{object} contains the fitted final
##'   learners when available.
##' @seealso \code{\link{run_rtmle}}, \code{\link{parse_learners}},
##'   \code{\link{learn_glm}}, \code{\link{learn_glmnet}},
##'   \code{\link{learn_ranger}}, \code{\link{learn_xgboost}}
##' @examples
##' # Note that the function is designed to be called from run_rtmle.
##' if (requireNamespace("ranger", quietly = TRUE)) {
##' data(simulated_cohort)
##' ld <- register_format(simulated_cohort)
##' A0 <- ld$timevar_data$A[date == 0, list(id, A_0 = value)]
##' d <- merge(ld$baseline_data, A0, by = "id")
##' fit <- superlearn(folds = 2,seed = 8,
##'                   learners = list("learn_glm",
##'                                   "my_ranger" = list(fun = "learn_ranger",
##'                                                      num.trees = 5),
##'                                   "learn_glmnet" = list(alpha = 1,
##'                                                      fun="learn_glmnet")),
##'                   character_formula = "A_0~sex+age+SBP",
##'                   outcome_variable = "A_0",
##'                   id_variable = "id",
##'                   data = d,
##'                   intervened_data = d)
##' }
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
superlearn <- function(folds,
                       seed,
                       learners,
                       parse_learners = TRUE,
                       character_formula,
                       outcome_variable,
                       outcome_variable_name,
                       id_variable,
                       data,
                       intervened_data,
                       ensemble_method = "discrete",
                       diagnostics,
                       ...){
    expected_levels = NULL
    N <- NROW(data)
    if (parse_learners){
        stopifnot(length(learners)>1)
        learners <- parse_learners(list(folds = folds,learners = learners))$learners
    }
    # set the seed
    if (!missing(seed) && is.numeric(seed)) set.seed(seed)
    # split the data
    random_split <- sample(1:folds,size = N,replace = TRUE)
    # count the levels of all factor variables and check for levels that do not occur in all folds
    if (length(char_variables <- names(Filter(is.character, data)))>0){
        stop(paste0("For superlearning all character predictor variables need to be converted to factor.\n",
                    "Offending variables: ",paste0(char_variables,collapse = ", ")))}
    factor_cols <- names(Filter(is.factor, data))
    if (length(factor_cols)>0){
        # need only pay attention to variables with more than 2 levels
        # because 2 level variables are caught by constant variable checks
        # when 1 level is missing
        multi_factor_levels <- data[, list(factor = factor_cols,
                                           expected_levels = sapply(.SD,uniqueN)),
                                    .SDcols = factor_cols][expected_levels>2]
    }else{
        multi_factor_levels <- NULL
    }
    # create empty level-1 data
    level_one_data <- data.table::data.table(id = 1:N,data[[outcome_variable]])
    setnames(level_one_data,c(id_variable,outcome_variable))
    # initialize a column in the level-1 data
    learner_names <- sapply(learners,function(x){x$name})
    for (this_learner_name in learner_names){
        set(level_one_data,j = this_learner_name,value = numeric(N))
    }
    # loop across folds
    for (k in 1:folds){
        character_formula_k <- character_formula
        # current learning data
        learn_data_k <- data[random_split != k]
        # current test data where treatment variables have values set under the intervention
        test_data_k <- data[random_split == k]
        # check for constant outcome and predictor variables
        current_constants <- sapply(learn_data_k, function(x){length(na.omit(unique(x)))==1})
        current_constants <- names(current_constants[current_constants])
        # CASE: no variation in the outcome variable
        #       predict the unique outcome value
        if (outcome_variable %in% current_constants){
            outcome_value_k <- na.omit(unique(learn_data_k[[outcome_variable]]))
            predicted_k <- rep(outcome_value_k,NROW(test_data_k))
            for (this_learner_name in learner_names) {
                set(level_one_data,j = this_learner_name,i = which(random_split == k),value = predicted_k)
            }
        } else {
            # CASE: no variation in one of the predictor variables
            #       remove these
            learn_data_k <- learn_data_k[,!(names(learn_data_k)%in%current_constants),with = FALSE]
            # now check if all multi_factor levels have all levels in learning set
            # fine if they do not have all levels in intervened_data
            if (NROW(multi_factor_levels)>0){
                vars_missing_levels_k <- unique(rbind(
                    multi_factor_levels,
                    data[, list(factor = multi_factor_levels$factor,
                                expected_levels = sapply(.SD,uniqueN)),
                         .SDcols = multi_factor_levels$factor]))[duplicated(factor)]$factor
                if (length(vars_missing_levels_k)>0){
                    if (!missing(diagnostics)){
                        diagnostics$missing_levels <- rbind(diagnostics$missing_levels,
                                                            data.table(Function = 'rtmle::superlearn',
                                                                       Outcome = outcome_variable_name,
                                                                       Event = "Not all levels in all folds",
                                                                       Info = paste0(vars_missing_levels_k,collapse = ", ")))
                    }
                    character_formula_k <- delete_variables_from_formula(character_formula = character_formula_k,
                                                                         delete_vars = vars_missing_levels_k)
                }
            }
            # remove constant predictor variables
            if (length(current_constants)>0){
                if (!missing(diagnostics)){
                    diagnostics$constant_predictors <- rbind(diagnostics$constant_predictors,
                                                             data.table(Function = 'rtmle::superlearn',
                                                                        Outcome = outcome_variable_name,
                                                                        Event = "constant_predictors",
                                                                        Fold = k,
                                                                        Info = paste0(current_constants,collapse = ", ")))
                }
                character_formula_k <- delete_variables_from_formula(character_formula = character_formula_k,
                                                                     delete_vars = current_constants)
                number_rhs_variables_k <- attr(character_formula_k,"number_rhs_variables")
            }else{
                number_rhs_variables_k <- length(attr(stats::terms(stats::formula(character_formula_k)),"term.labels"))
            }
            # if there are no covariates then we simply predict the mean 
            if (number_rhs_variables_k == 0){
                # here we assume that the outcome is binary or a predicted continuous value 
                predicted_k <- rep(mean(learn_data_k[[outcome_variable]],na.rm = TRUE),NROW(test_data_k))
                for (this_learner_name in learner_names) {
                    set(level_one_data,
                        j = this_learner_name,
                        i = which(random_split == k),
                        value = as.vector(predicted_k))
                }
            }else{
                for (this_learner in learners) {
                    learner_args <- c(list(character_formula = character_formula_k,
                                           data = learn_data_k,
                                           intervened_data = test_data_k),
                                      ## learner specific arguments such as tuning parameters
                                      this_learner$args)
                    if (inherits(try( 
                        predicted_k <- as_learner_output(do.call(this_learner$fun, learner_args))
                    ),"try-error")){
                        stop(paste0("Learning failed in fold ",k," with learner ",this_learner))
                    }
                    set(level_one_data,
                        j = this_learner$name,
                        i = which(random_split == k),
                        value = as.vector(predicted_k$predicted_values))
                }
            }
        }
    }
    ensemble_weights <- switch(tolower(ensemble_method),
                               # Non-Negative Least Squares
                               "nnls" = {
                                   L1 = as.matrix(level_one_data[,learner_names,with = FALSE])
                                   nnls_solution <- nnls::nnls(A = L1,b = level_one_data[[outcome_variable]])$x
                                   nnls_weights <- nnls_solution / sum(nnls_solution)
                                   names(nnls_weights) <- learner_names
                                   nnls_weights
                               },
                               "ipa" = {
                                   null_model_Brier_score <- mean((level_one_data[[outcome_variable]]-mean(data[[outcome_variable]]))^2)
                                   Brier_score <- sapply(learner_names,function(this_learner_name){
                                       mean(
                                       (level_one_data[[outcome_variable]]-level_one_data[[this_learner_name]]
                                       )^2,na.rm = TRUE)})
                                   ipa_weights <- 1-(Brier_score/null_model_Brier_score)
                                   # remove models that perform worse than the null model
                                   ipa_weights <- pmax(0,ipa_weights)
                                   if (any(ipa_weights>0)){
                                       ipa_weights <- ipa_weights/sum(ipa_weights)
                                   }
                                   names(ipa_weights) <- learner_names
                                   ipa_weights                                       
                               }, {
                                   # default is discrete super learner (for now)
                                   # choose the minimizer of the Brier score
                                   x <- sapply(learner_names,function(this_learner_name){
                                       mean(
                                       (level_one_data[[outcome_variable]]-level_one_data[[this_learner_name]]
                                       )^2,na.rm = TRUE)})
                                   names(x) <- learner_names
                                   x <- sort(x)
                                   discrete_weights <- c(1,rep(0,length(x)-1))
                                   names(discrete_weights) <- names(x)
                                   discrete_weights
                               })
    ## winner_name <- names(x)[1]
    ## names(learners) <- learner_names
    ## winner <- learners[[winner_name]]
    names(learners) <- learner_names
    learner_args <- list(character_formula = character_formula,
                         data = data,
                         intervened_data = intervened_data)
    if (all(ensemble_weights == 0)){
        # NULL model prediction
        if (!missing(diagnostics)){
            diagnostics$superlearner_null_model <- rbind(diagnostics$superlearner_null_model,
                                                         data.table(Function ='rtmle::superlearn',
                                                                    Outcome = outcome_variable_name,
                                                                    Event = "All learners performed worse than the null model"))
        }
        predicted_values <- rep(mean(data[[outcome_variable]]),NROW(intervened_data))
        fitted_objects <- NULL
    }else{
        predicted_value_ensemble <- lapply(names(ensemble_weights)[ensemble_weights>0],function(this_learner){
            ## add winner's arguments but do not include duplicate arguments
            predicted_values_args <- c(learners[[this_learner]]$args,
                                       learner_args[!names(learner_args)%in%names(learners[[this_learner]]$args)])
            as_learner_output(do.call(what = learners[[this_learner]]$fun,predicted_values_args))
        })
        if (length(predicted_value_ensemble) == 1){
            predicted_values <- predicted_value_ensemble[[1]]$predicted_values
        }else{
            predicted_values <- do.call(cbind,lapply(predicted_value_ensemble,`[[`,"predicted_values"))%*%ensemble_weights[ensemble_weights>0]
        }
        fitted_objects <- lapply(predicted_value_ensemble,`[[`,"object")
        names(fitted_objects) <- names(ensemble_weights)[ensemble_weights>0]
    }
    learner_output(predicted_values = predicted_values,
                   diagnostics = if (!missing(diagnostics)) diagnostics else NULL,
                   fit = ensemble_weights,
                   object = list(ensemble_weights = ensemble_weights,
                                 fitted_objects = fitted_objects))
}


######################################################################
### superlearn.R ends here
