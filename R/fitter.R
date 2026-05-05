### fitter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: jan 25 2026 (09:26) 
## Version: 
## Last-Updated: maj  4 2026 (12:57) 
##           By: Thomas Alexander Gerds
##     Update #: 101
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
fitter <- function(intervention_node,
                   learner,
                   formula,
                   data,
                   intervened_data,
                   id_variable,
                   minority_threshold,
                   seed,
                   diagnostics = NULL,
                   save_fitted_objects = FALSE){
    current_data <- data.table::as.data.table(model.frame(formula = stats::formula(formula),data = data))
    ## I(1*(A==1)) has class 'AsIs' but the values are numeric, that's why as.numeric is applied
    current_outcome_name <- names(current_data)[[1]]
    outcome_variable <- as.numeric(current_data[[1]])
    names(current_data)[1] <- "rtmle_outcome"
    ff <- change_outcome_in_formula(formula,"rtmle_outcome")
    current_constants <- sapply(current_data, function(x){length(na.omit(unique(x)))<=1})
    n_outcome_values <- length(unique(current_data[[1]]))
    if (current_constants[[1]] == TRUE) {minority_count <- 0}
    if (current_constants[[1]] == TRUE ||
        # when the outcome is binary but has less than threshold many occurences of the minority group
        # simply predict the mean
        (n_outcome_values == 2 && (minority_count <- min(table(current_data[[1]])))<minority_threshold)){
        if (!missing(diagnostics)){
            diag_row <- data.table(Function ='rtmle::fitter',
                                   Intervention_node = intervention_node,
                                   Outcome = gsub("\"","'",current_outcome_name),
                                   Event = "Constant outcome or below minority_threshold")
            if (is.null(diagnostics)){
                diagnostics <- list(constant_outcome_variables = diag_row)
            }else{
                diagnostics$constant_outcome_variables <- rbind(diagnostics$constant_outcome_variables,diag_row)
            }
        }
        mean_Y <- mean(outcome_variable)
        predicted_values <- rep(mean_Y,NROW(intervened_data))
        intercept <- log(mean_Y/(1-mean_Y))
        predicted_values <- list(predicted_values = predicted_values,
                                 diagnostics = diagnostics,
                                 fit = matrix(intercept,ncol = 1,nrow = 1,dimnames = list("(Intercept)","Estimate")))
    }else{
        # remove all currently constant predictor variables
        # from the formula
        if (any(current_constants)) {
            current_constants <- names(current_constants[current_constants])
            current_data <- current_data[,!(names(current_data)%in%current_constants),with = FALSE]
            if (!missing(diagnostics)){
                diag_row <- data.table(
                    Function ='rtmle::fitter', 
                    Intervention_node = intervention_node,
                    Outcome = gsub("\"","'",current_outcome_name),
                    Event = "Constant predictors",
                    Info = paste0(current_constants,collapse = ", "))
                if (is.null(diagnostics)){
                    diagnostics <- list(constant_predictor_variables = diag_row) 
                }else{
                    diagnostics$constant_predictor_variables <- rbind(diagnostics$constant_predictor_variables,
                                                                      diag_row)
                }
            }
            ff <- delete_variables_from_formula(character_formula = ff,delete_vars = current_constants)
            number_rhs_variables <- attr(ff,"number_rhs_variables",exact = TRUE)
        }else{
            current_constants <- NULL
            number_rhs_variables <- length(attr(stats::terms(stats::formula(ff)),"term.labels",exact = TRUE))
        }
        if (number_rhs_variables == 0){
            mean_Y <- mean(outcome_variable)
            predicted_values <- rep(mean_Y,NROW(intervened_data))
            if (!missing(diagnostics)){
                diag_row <- data.table(
                    Function ='rtmle::fitter', 
                    Intervention_node = intervention_node,
                    Outcome = gsub("\"","'",current_outcome_name),
                    Event = "Empty rhs of formula")
                if (is.null(diagnostics)){
                    diagnostics <- list(empty_rhs = diag_row)
                }else{
                    diagnostics$empty_rhs <- rbind(diagnostics$empty_rhs,diag_row)
                }
            }
            intercept <- log(mean_Y/(1-mean_Y))
            predicted_values <- list(predicted_values = predicted_values,
                                     diagnostics = diagnostics,
                                     fit = matrix(intercept,ncol = 1,nrow = 1,dimnames = list("(Intercept)","Estimate")))
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
                               outcome_variable = "rtmle_outcome",
                               outcome_variable_name = current_outcome_name,
                               id_variable = id_variable))
                if (inherits(try(
                    predicted_values <- do.call("superlearn",
                                                c(args,
                                                  list(diagnostics = diagnostics,
                                                       current_constants,
                                                       seed = seed,
                                                       save_fitted_objects = save_fitted_objects))),
                    silent = FALSE),
                    "try-error")) {
                    stop(paste0("Failed to superlearn/crossfit with formula ",ff,"\nwhere the outcome is: ",
                                ifelse(current_outcome_name == "rtmle_predicted_outcome",
                                       paste0("the predicted outcome at intervention node (time): ",intervention_node),
                                       current_outcome_name)))
                }
            }else{
                #
                # this is a single learner 
                #
                learner_warnings <- list()
                result <- withCallingHandlers({
                    if (inherits(
                        try(
                            predicted_values <- do.call(
                                learner$fun,
                                c(args,c(learner$args,list(save_fitted_objects = save_fitted_objects)))
                            ),
                            silent = FALSE
                        ),
                        "try-error"
                    )) {
                        stop(paste0(
                            "Failed to superlearn/crossfit with formula ", ff,
                            "\nwhere the outcome is: ",
                            ifelse(
                                current_outcome_name == "rtmle_predicted_outcome",
                                paste0("the predicted outcome at intervention node (time): ", intervention_node),
                                current_outcome_name
                            )
                        ))
                    }
                    predicted_values
                },
                warning = function(w) {
                    learner_warnings <<- c(learner_warnings, conditionMessage(w))
                    invokeRestart("muffleWarning")
                })
                if (length(learner_warnings)>0){
                    if (length(learner_warnings)>1){
                        learner_warnings <- paste0("warning 1: ",learner_warnings[1]," warning 2: ",learner_warnings[2]," ...")
                    }
                    predicted_values$diagnostics <- c(predicted_values$diagnostics,
                                                      list(learner_warnings = data.table(learner = learner$name,
                                                                                         warning = learner_warnings)))
                }
            }
        }
    }
    predicted_values
}

######################################################################
### fitter.R ends here
