### superlearn.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 31 2024 (07:29) 
## Version: 
## Last-Updated: Apr 11 2025 (16:40) 
##           By: Thomas Alexander Gerds
##     Update #: 114
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
                       outcome_target_level,
                       id_variable,
                       data,
                       intervened_data,
                       ensemble,
                       ...){
    N <- NROW(data)
    stopifnot(length(learners)>1)
    learners <- parse_learners(learners)
    if (!missing(seed) && is.numeric(seed)) set.seed(seed)
    split <- sample(1:folds,size = N,replace = TRUE)
    level_one_data <- data.table::data.table(id = 1:N,data[[outcome_variable]])
    setnames(level_one_data,c(id_variable,outcome_variable))
    # initialize a column in the level_one_data
    for (this_learner_name in names(learners))
        set(level_one_data,j = this_learner_name,value = numeric(N))
    for (k in 1:folds){
        learn_data <- data[split == k]
        i_data <- intervened_data[split != k]
        learning_args <- list(character_formula = character_formula,
                              data = learn_data,
                              intervened_data = i_data)
        # when there is no variation in the outcome variable then the predicted risk of all learners
        # is set to this unique value of the outcome variable
        if (length(unival <- unique(na.omit(learn_data[[outcome_variable]]))) == 1){
            if (!is.null(outcome_target_level)){
                predicted_k <- rep(1*(outcome_target_level == unival),sum(split != k))
            }else {
                predicted_k <- rep(unival,sum(split != k))
            }
            for (this_learner_name in names(learners)) {
                set(level_one_data,
                    j = this_learner_name,
                    i = which(split != k),
                    value = predicted_k)
            }
        } else{
            for (this_learner_name in names(learners)) {
                this_learner <- learners[[this_learner_name]]
                learner_args <- list(character_formula = character_formula,
                                     data = learn_data,
                                     intervened_data = i_data)
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
                    i = which(split != k),
                    value = as.vector(predicted_k))
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
