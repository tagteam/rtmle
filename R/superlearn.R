### superlearn.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 31 2024 (07:29) 
## Version: 
## Last-Updated: Nov  1 2024 (07:29) 
##           By: Thomas Alexander Gerds
##     Update #: 38
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
    for (l in learners) set(level_one_data,j = l,value = numeric(N))
    for (k in 1:folds){
        ## print(k)
        learn_data <- data[split == k]
        i_data <- intervened_data[split != k]
        for (l in learners) {
            ## print(l)
            if (length(unique(learn_data[[outcome_variable]])) == 1){
                if (length(grep("Censored_",outcome_variable)>0)) {
                    predicted_k <- rep(1*("uncensored" == unique(learn_data[[outcome_variable]])),sum(split != k))
                }else {
                    predicted_k <- rep(unique(learn_data[[outcome_variable]]),sum(split != k))
                }
            } else
                predicted_k <- do.call(l,list(character_formula = character_formula,data = learn_data,intervened_data = i_data,...))
            set(level_one_data,
                j = l,
                i = which(split != k),
                value = predicted_k)
        }
    }
    ## print(outcome_variable)
    ## if (outcome_variable == "primary.outcome_3") browser(skipCalls=1L)
    # discrete super learner (for now)
    if (length(grep("Censored_",outcome_variable)>0)){
        x <- sapply(learners,function(l){mean(((1*(level_one_data[[outcome_variable]] == "uncensored"))-level_one_data[[l]])^2)})
    } else{
        # FIXME: NA values of the outcome should not make it until here
        x <- sapply(learners,function(l){mean((level_one_data[[outcome_variable]]-level_one_data[[l]])^2,na.rm = TRUE)})
    }
    winner <- names(x)[which.min(x)]
    ## print(winner)
    return(level_one_data[[winner]])
}


######################################################################
### superlearn.R ends here
