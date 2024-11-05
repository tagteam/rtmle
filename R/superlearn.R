### superlearn.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 31 2024 (07:29) 
## Version: 
## Last-Updated: Nov  5 2024 (15:27) 
##           By: Thomas Alexander Gerds
##     Update #: 57
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
superlearn <- function(folds,
                       seed,
                       learners,
                       character_formula,
                       outcome_variable,
                       id_variable,
                       data,
                       intervened_data,
                       ...){
    N <- NROW(data)
    if (!missing(seed) && is.numeric(seed)) set.seed(seed)
    split <- sample(1:folds,size = N,replace = TRUE)
    level_one_data <- data.table::data.table(id = 1:N,data[[outcome_variable]])
    setnames(level_one_data,c(id_variable,outcome_variable))
    for (this_learner_name in names(learners)) set(level_one_data,j = this_learner_name,value = numeric(N))
    for (k in 1:folds){
        ## print(k)
        learn_data <- data[split == k]
        i_data <- intervened_data[split != k]
        learning_args <- list(character_formula = character_formula,
                              data = learn_data,
                              intervened_data = i_data)
        for (this_learner_name in names(learners)) {
            this_learner <- learners[[this_learner_name]]
            if (length(unique(na.omit(learn_data[[outcome_variable]]))) == 1){
                if (length(grep("Censored_",outcome_variable)>0)) {
                    predicted_k <- rep(1*("uncensored" == unique(learn_data[[outcome_variable]])),sum(split != k))
                }else {
                    predicted_k <- rep(unique(learn_data[[outcome_variable]]),sum(split != k))
                }
            } else{
                learner_args <- list(character_formula = character_formula,
                                     data = learn_data,
                                     intervened_data = i_data)
                ## adding learner specific tuning parameters
                learner_fun_pos <- match("learner_fun",names(learners[[this_learner_name]]))
                if (length(this_learner) == 1)
                    this_learner_args <- NULL
                else
                    this_learner_args <- this_learner[[-learner_fun_pos]]
                if (is.na(learner_fun_pos)) stop(paste0("Cannot find learner function for learner ",this_learner_name))
                if (inherits(try( 
                    predicted_k <- do.call(this_learner[[learner_fun_pos]],
                                           c(learner_args,this_learner_args))
                ),"try-error")){
                    ## browser(skipCalls=1L)
                    stop(paste0("Learning failed in fold ",k," with learner ",this_learner))
                }
            }
            set(level_one_data,
                j = this_learner_name,
                i = which(split != k),
                value = predicted_k)
        }
    }
    ## print(outcome_variable)
    ## if (outcome_variable == "primary.outcome_3") browser(skipCalls=1L)
    # discrete super learner (for now)
    if (length(grep("Censored_",outcome_variable)>0)){
        x <- sapply(names(learners),function(this_learner_name){mean(((1*(level_one_data[[outcome_variable]] == "uncensored"))-level_one_data[[this_learner_name]])^2)})
    } else{
        # FIXME: NA values of the outcome should not make it until here
        x <- sapply(names(learners),function(this_learner_name){mean((level_one_data[[outcome_variable]]-level_one_data[[this_learner_name]])^2,na.rm = TRUE)})
    }
    winner <- names(x)[which.min(x)]
    print(winner)
    return(level_one_data[[winner]])
}


######################################################################
### superlearn.R ends here
