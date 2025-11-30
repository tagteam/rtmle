### superlearn.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 31 2024 (07:29) 
## Version: 
## Last-Updated: nov 30 2025 (13:16) 
##           By: Thomas Alexander Gerds
##     Update #: 159
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Super learning nuisance parameter models for longitudinal TMLE analyses
##'
##' The function collaborates with \code{run_rtmle}. It 
##' It uses cross-fitting to produce level-one data, where the predicted probability of each learner
##' is a column which can be compared with the outcome. The discrete super learner chooses the
##' learner lowest quadratic prediction error. The ensemble super learner combines the learners.
##' 
##' @title Super learning for rtmle
##' @param folds Number of folds
##' @param seed Random seed
##' @param learners List of learners. Each learner can either be a character string or a list. If the learner is a
##' character string it must be the name of a learner function, such as \code{\link{learn_glm}},
##'  \code{\link{learn_glmnet}},  \code{\link{learn_ranger}}. If the learner is a list then either the name of the learner
##' is a learner function or the list has an element \code{learner_fun} which is the name of a learner function. The other elements of the list are
##' passed on as additional arguments to the learner function.
##' @param character_formula A formula (passed as a character string!) to parse the outcome and the predictor variables. 
##' @param outcome_variable The name of the outcome variable (same as \code{all.vars(formula(character_formula))[[1]]} provided to avoid overhead in the parsing of the formula.
##' @param outcome_target_level The level of the binary outcome variable for which the superlearner predicts the risk.  
##' @param id_variable The name of the subject id variable.
##' @param data The data for learning.
##' @param intervened_data The data were all intervention variables are readily set according to the intervention protocol.
##' @param ensemble Type of superlearning. When \code{"winner"} the discrete super learner is the learner with the lowest
##' level-one mean squared prediction error (Brier score).
##' @param ... Not (yet) used. 
##' @return The level-one data column of the discrete superlearner. 
##' @seealso \code{\link{run_rtmle}}, \code{learn_glm}, \code{learn_glmnet}, \code{learn_ranger}
##' @examples
##' # Note that the function is designed to be called from run_rtmle.  
##' library(ranger)
##' library(glmnet)
##' set.seed(17)
##' ld <- simulate_long_data(n = 113,
##'                          number_visits = 20,
##'                          beta = list(A_on_Y = -.2,
##'                                      A0_on_Y = -0.3,A0_on_A = 6),
##'                          register_format = FALSE)
##' fit <- superlearn(folds = 2,seed = 8,
##'                   learners = list("learn_glm",
##'                                   "my_ranger" = list(learner_fun = "learn_ranger",
##'                                                      num.trees = 5),
##'                                   "learn_glmnet" = list(alpha = 1)),
##'                   character_formula = "A_0~L_0+sex+age",
##'                   outcome_variable = "A_0",
##'                   outcome_target_level="1",
##'                   id_variable = "id",
##'                   data = ld[!duplicated(id)],
##'                   intervened_data = ld[!duplicated(id)])
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
superlearn <- function(folds,
                       seed,
                       learners,
                       character_formula,
                       outcome_variable,
                       outcome_target_level, # for treatment and censoring variables
                       id_variable,
                       data,
                       intervened_data,
                       ensemble,
                       ...){
    N <- NROW(data)
    stopifnot(length(learners)>1)
    learners <- parse_learners(learners)
    # set the seed
    if (!missing(seed) && is.numeric(seed)) set.seed(seed)
    # split the data
    random_split <- sample(1:folds,size = N,replace = TRUE)
    # count the levels of all factor variables and check for levels that do not occur in all folds
    if (length(char_variables <- names(Filter(is.character, data)))>0){
        stop(paste0("For save machine learning all character predictor variables need to be converted to factor.\n",
                    "Offending variables: ",paste0(char_variables,collapse = ", ")))}
    factor_cols <- names(Filter(is.factor, data))
    if (length(factor_cols)>0){
        # need only pay attention to variables with more than 2 levels
        # because 2 level variables are caught by constant variable checks
        # when 1 level is missing
        multi_factor_levels <- data[, .(factor = factor_cols,
                                        expected_levels = sapply(.SD,uniqueN)),
                                    .SDcols = factor_cols][expected_levels>2]
    }
    # create emtpy level-1 data
    level_one_data <- data.table::data.table(id = 1:N,data[[outcome_variable]])
    setnames(level_one_data,c(id_variable,outcome_variable))
    # initialize a column in the level-1 data
    for (this_learner_name in names(learners)){
        set(level_one_data,j = this_learner_name,value = numeric(N))
    }
    # loop across folds
    for (k in 1:folds){
        character_formula_k <- character_formula
        # current learning data
        learn_data_k <- data[random_split != k]
        # current test data where treatment variables have values set under the intervention
        i_data_k <- intervened_data[random_split == k]
        # check for constant outcome and predictor variables
        current_constants <- sapply(learn_data_k, function(x){length(na.omit(unique(x)))==1})
        current_constants <- names(current_constants[current_constants])
        # CASE: no variation in one of the treatment variables
        # FIXME: should we report this or just stop the algorithm?
        # CASE: no variation in the outcome variable
        #       predict the unique outcome value
        if (outcome_variable %in% current_constants){
            if (!is.null(outcome_target_level)){
                outcome_value_k <- na.omit(unique(learn_data_k[[outcome_variable]]))
                predicted_k <- rep(1*(outcome_target_level == outcome_value_k),NROW(i_data_k))
            }else {
                predicted_k <- rep(outcome_value_k,NROW(i_data_k))
            }
            for (this_learner_name in names(learners)) {
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
                    data[, .(factor = multi_factor_levels$factor,
                             expected_levels = sapply(.SD,uniqueN)),
                         .SDcols = multi_factor_levels$factor]))[duplicated(factor)]$factor
                if (length(vars_missing_levels_k)>0){
                    warning(paste0("Outcome: ",outcome_variable,": the following variables have not all levels in all folds: ",
                                   paste0(vars_missing_levels_k,collapse = ", ")))
                    character_formula_k <- delete_variables_from_formula(character_formula = character_formula_k,
                                                                         delete_vars = vars_missing_levels_k)
                }
            }


            # remove constant predictor variables
            if (length(current_constants)>0){
                warning(paste0("Outcome: ",outcome_variable,": the following variables are constant in fold ",k,": ",
                               paste0(current_constants,collapse = ", ")))
                character_formula_k <- delete_variables_from_formula(character_formula = character_formula_k,
                                                                     delete_vars = current_constants)
                number_rhs_variables_k <- attr(character_formula_k,"number_rhs_variables")
            }else{
                number_rhs_variables_k <- length(attr(stats::terms(stats::formula(character_formula_k)),"term.labels"))
            }
            # if there are no covariates then we simply predict the mean 
            if (number_rhs_variables_k == 0){
                # here we assume that the outcome is binary or a predicted continuous value 
                if (!is.null(outcome_target_level)){
                    predicted_k <- rep(mean(outcome_target_level == learn_data_k[[outcome_variable]],na.rm = TRUE),NROW(i_data_k))
                }else {
                    predicted_k <- rep(mean(learn_data_k[[outcome_variable]],na.rm = TRUE),NROW(i_data_k))
                }
                for (this_learner_name in names(learners)) {
                    set(level_one_data,
                        j = this_learner_name,
                        i = which(random_split == k),
                        value = as.vector(predicted_k))
                }
            }else{
                for (this_learner_name in names(learners)) {
                    this_learner <- learners[[this_learner_name]]
                    learner_args <- list(character_formula = character_formula_k,
                                         data = learn_data_k,
                                         intervened_data = i_data_k)
                    ## adding learner specific tuning parameters
                    learner_fun_pos <- match("learner_fun",
                                             names(learners[[this_learner_name]]))
                    if (length(this_learner) == 1){
                        this_learner_args <- NULL
                    } else{
                        ## if (is.na(learner_fun_pos)) stop(paste0("Cannot find learner function for learner ",this_learner_name))
                        this_learner_args <- this_learner[-learner_fun_pos]
                    }
                    if (inherits(try( 
                        predicted_k <- do.call(this_learner[[learner_fun_pos]],
                                               c(learner_args,this_learner_args))
                    ),"try-error")){
                        stop(paste0("Learning failed in fold ",k," with learner ",this_learner))
                    }
                    set(level_one_data,
                        j = this_learner_name,
                        i = which(random_split == k),
                        value = as.vector(predicted_k))
                }
            }
        }
    }
    # discrete super learner (for now)
    # choose the minimizer of the Brier score
    if (!is.null(outcome_target_level)){
        x <- sapply(names(learners), function(this_learner_name){
            mean(((1*(level_one_data[[outcome_variable]] == outcome_target_level))-level_one_data[[this_learner_name]])^2)
        })
    } else{
        # FIXME: NA values of the outcome should not make it until here
        x <- sapply(names(learners),function(this_learner_name){
            outcome <- level_one_data[[outcome_variable]]
            mean((outcome-level_one_data[[this_learner_name]])^2,na.rm = TRUE)})
    }
    winner <- names(x)[which.min(x)]
    ## print(winner)
    return(level_one_data[[winner]])
}


######################################################################
### superlearn.R ends here
