### fitter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: jan 25 2026 (09:26) 
## Version: 
## Last-Updated: mar 25 2026 (17:34) 
##           By: Thomas Alexander Gerds
##     Update #: 82
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
                   diagnostics = NULL){
    current_data <- data.table::as.data.table(model.frame(formula = stats::formula(formula),data = data))
    ## I(1*(A==1)) has class 'AsIs' but the values are numeric, that's why as.numeric is applied
    current_outcome_name <- names(current_data)[[1]]
    outcome_variable <- as.numeric(current_data[[1]])
    names(current_data)[1] <- "rtmle_outcome"
    ff <- change_outcome_in_formula(formula,"rtmle_outcome")
    current_constants <- sapply(current_data, function(x){length(na.omit(unique(x)))<=1})
    n_outcome_values <- length(unique(current_data[[1]]))
    if (current_constants[[1]] == TRUE ||
        # when the outcome is binary but has less than threshold many occurences of the minority group
        # simply predict the mean
        (n_outcome_values == 2 && any(table(current_data[[1]])<minority_threshold))){
        if (!missing(diagnostics)){
            if (is.null(diagnostics)){
                diagnostics <- list(constant_outcome_variables <- c(
                                        diagnostics$constant_outcome_variables,
                                        paste0("Intervention node ",intervention_node,", outcome ",current_outcome_name,": is constant or the variation is below 'x$tuning_parameters$minority_threshold'.")
                                    ))
            } else{
                diagnostics$constant_outcome_variables <- c(
                    diagnostics$constant_outcome_variables,
                    paste0("Intervention node ",intervention_node,", outcome ",current_outcome_name,": is constant or the variation is below 'x$tuning_parameters$minority_threshold'.")
                )
            }
        }
        mean_Y <- mean(outcome_variable)
        predicted_values <- rep(mean_Y,NROW(intervened_data))
        intercept <- log(mean_Y/(1-mean_Y))
        data.table::setattr(predicted_values,"fit",matrix(intercept,ncol = 1,nrow = 1,dimnames = list("(Intercept)","Estimate")))
        data.table::setattr(predicted_values,"diagnostics",diagnostics)
    }else{
        # remove all currently constant predictor variables
        # from the formula
        if (any(current_constants)) {
            current_constants <- names(current_constants[current_constants])
            current_data <- current_data[,!(names(current_data)%in%current_constants),with = FALSE]
            if (!missing(diagnostics)){
                if (is.null(diagnostics)){
                    diagnostics <- list(constant_predictor_variables = paste0("Intervention node ",intervention_node,", outcome ",current_outcome_name,": the following predictor variables are constant: ",
                                                                              paste0(current_constants,collapse = ", ")))
                }else{
                    diagnostics$constant_predictor_variables <- c(diagnostics$constant_predictor_variables,
                                                                  paste0("Intervention node ",intervention_node,", outcome ",current_outcome_name,": the following predictor variables are constant: ",
                                                                         paste0(current_constants,collapse = ", ")))
                }
            }
            ff <- delete_variables_from_formula(character_formula = ff,delete_vars = current_constants)
            number_rhs_variables <- attr(ff,"number_rhs_variables")
        }else{
            current_constants <- NULL
            number_rhs_variables <- length(attr(stats::terms(stats::formula(ff)),"term.labels"))
        }
        if (number_rhs_variables == 0){
            mean_Y <- mean(outcome_variable)
            predicted_values <- rep(mean_Y,NROW(intervened_data))
            if (!missing(diagnostics)){
                if (is.null(diagnostics)){
                    diagnostics <- list(empty_rhs = paste0("Intervention node ",intervention_node,", outcome ",current_outcome_name))
                }else{
                    diagnostics$empty_rhs <- c(diagnostics$empty_rhs,
                                               paste0("Intervention node ",intervention_node,", outcome ",current_outcome_name))
                }
                data.table::setattr(predicted_values,"diagnostics",diagnostics)
            }
            intercept <- log(mean_Y/(1-mean_Y))
            data.table::setattr(predicted_values,"fit",matrix(intercept,ncol = 1,nrow = 1,dimnames = list("(Intercept)","Estimate")))
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
                                                       seed = seed))),silent = FALSE),
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
                if (inherits(try(
                    predicted_values <- do.call(learner$fun,c(args,learner$args)),silent = FALSE),
                    "try-error")) {
                    stop(paste0("Failed to superlearn/crossfit with formula ",ff,"\nwhere the outcome is: ",
                                ifelse(current_outcome_name == "rtmle_predicted_outcome",
                                       paste0("the predicted outcome at intervention node (time): ",intervention_node),
                                       current_outcome_name)))
                }else{
                    if (!missing(diagnostics) && !is.null(diagnostics)){
                        if (is.null(attr(predicted_values,"diagnostics"))){
                            data.table::setattr(predicted_values,"diagnostics",diagnostics)
                        }
                    }
                }
            }
        }
    }
    predicted_values
}

######################################################################
### fitter.R ends here
