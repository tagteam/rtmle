### fitter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: jan 25 2026 (09:26) 
## Version: 
## Last-Updated: jan 30 2026 (11:47) 
##           By: Thomas Alexander Gerds
##     Update #: 21
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
fitter <- function(learner,
                   formula,
                   data,
                   intervened_data,
                   id_variable){
    ## if (any(is.na(data))){
    ## has_missing <- sapply(data,function(x)sum(is.na(x)))
    ## has_missing <- has_missing[has_missing != 0]
    ## stop(paste0("Missing values detected in data for fitting nuisance parameter models at time ",k,":\n",
    ## paste0(names(has_missing),": n=",has_missing,collapse = "\n")))
    ## }
    current_data <- model.frame(formula = stats::formula(formula),data = data)
    outcome_variable <- current_data[[1]]
    names(current_data)[1] <- "rtmle_binary"
    ff <- change_outcome_in_formula(formula,"rtmle_binary")
    current_constants <- sapply(current_data, function(x){length(na.omit(unique(x)))<=1})
    if (current_constants[[1]] == TRUE){
        predicted_values <- outcome_variable
        attr(predicted_values,"fit") <- paste0("No variation in this variable.")
    }else{
        # remove all currently constant predictor variables
        # from the formula
        if (any(current_constants)) {
            current_constants <- names(current_constants[current_constants])
            current_data <- current_data[,!(names(current_data)%in%current_constants),with = FALSE]
            ff <- delete_variables_from_formula(character_formula = ff,delete_vars = current_constants)
            number_rhs_variables <- attr(formula,"number_rhs_variables")
        }else{
            current_constants <- NULL
            number_rhs_variables <- length(attr(stats::terms(stats::formula(ff)),"term.labels"))
        }
        if (number_rhs_variables == 0){
            predicted_values <- rep(mean(outcome_variable),NROW(current_data))
            attr(predicted_values,"fit") <- paste0("No covariates. Predicted average value.")
        } else {
            args <- list(character_formula = ff,
                         data = current_data,
                         intervened_data = intervened_data[,names(current_data)[-1],with = FALSE])
            if (learner$name == "superlearn"){
                #
                # this needs a list of learners
                #
                args <- c(args,
                          learner$args,  # folds and other arguments for superlearning
                          list(learners = learner$learners,
                               parse_learners = FALSE,
                               outcome_variable = outcome_variable,
                               id_variable = id_variable))
                if (inherits(try(
                    predicted_values <- do.call("superlearn",c(args,list(seed = seed))),silent = FALSE),
                    "try-error")) {
                    stop(paste0("Failed to superlearn/crossfit with formula ",ff))
                }
            }else{
                #
                # this is a single learner 
                #
                if (inherits(try(
                    predicted_values <- do.call(learner$fun,c(args,learner$args)),silent = FALSE),
                    "try-error")) {
                    stop(paste0("Failed to learn/predict with formula ",ff))
                }
            }
        }
    }
    predicted_values
}

######################################################################
### fitter.R ends here
