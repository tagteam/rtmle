### learn_xgboost.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 28 2024 (09:26) 
## Version: 
## Last-Updated: apr 29 2026 (07:35) 
##           By: Thomas Alexander Gerds
##     Update #: 137
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Learn nuisance-parameter models with xgboost
##'
##' Learns nuisance-parameter models for TMLE and predicts probabilities in
##' intervention-updated data using \code{\link[xgboost]{xgboost}}.
##'
##' @title Nuisance-parameter learner based on xgboost
##' @param character_formula Formula for the nuisance parameter, supplied as a
##'   character string.
##' @param data Data used for learning.
##' @param intervened_data Data used for prediction after intervention variables
##'   have been set according to a protocol.
##' @param ... Additional arguments passed to \code{\link[xgboost]{xgboost}},
##'   including hyperparameters.
##' @return A vector of predicted probabilities with the fitted model stored in
##'   the \code{"fit"} attribute.
##' @seealso \code{\link{superlearn}}, \code{\link{learn_glm}},
##'   \code{\link{learn_glmnet}}, \code{\link{learn_ranger}}
##' @examples
#' tau <- 3
#' data(simulated_cohort)
#' ld <- register_format(simulated_cohort)
#' x <- rtmle_init(time_grid = seq(0,20,4),name_id = "id",
#'                 name_outcome = "stroke",name_competing = "death",
#'                 name_censoring = "dropout",censored_label = "censored")
#' x <- add_long_data(x,
#'                    outcome_data=ld$timevar_data$stroke[!duplicated(id)],
#'                    censored_data=ld$timevar_data$dropout,
#'                    competing_data=ld$timevar_data$death,
#'                    timevar_data=ld$timevar_data[c("bleeding","changeSBP","A","B")])
#' x <- add_baseline_data(x,data=ld$baseline_data)
#' x <- long_to_wide(x,start_followup_date=0)
#' x <- prepare_rtmle_data(x)
#' x <- protocol(x,name = "Always_A",
#'                     intervention = data.frame(time=x$intervention_nodes,
#'                                               "A" = factor("1",levels = c("0","1"))))
#' x <- protocol(x,name = "Never_A",
#'                     intervention = data.frame(time=x$intervention_nodes,
#'                                               "A" = factor("0",levels = c("0","1"))))
#' x <- target(x,name = "Outcome_risk",
#'                   estimator = "tmle",
#'                   protocols = c("Always_A","Never_A"))
#' x <- model_formula(x)
#' x <- run_rtmle(x,learner = "learn_xgboost",time_horizon = 3)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
learn_xgboost <- function(character_formula,
                          data,
                          intervened_data,
                          ...){
    if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("xgboost library required")
    }
    control <- prodlim::SmartControl(call= list(...),
                                     keys=c("xgboost"),
                                     ignore=NULL,
                                     ignore.case=TRUE,
                                     defaults=list("xgboost"=list(max_depth = 2L,
                                                                  learning_rate = 1.0,
                                                                  nrounds = 2L,
                                                                  subsample = 1.0,
                                                                  reg_lambda = 1,
                                                                  objective = "reg:squarederror")),
                                     forced=NULL,
                                     verbose=TRUE)
    # extract the data
    sf <- Publish::specialFrame(stats::formula(character_formula),
                                data = data,
                                specials = NULL,
                                na.action = na.omit)
    Y <- as.numeric(sf$response[[1]])
    ## d <- xgboost::xgb.DMatrix(sf$design, label = Y)
    if (!inherits(try(
             fit <- do.call(xgboost::xgboost,c(list(x = sf$design,y = Y),control$xgboost))
            ,silent = TRUE),
             "try-error")){
    }else{
        stop(paste0("\nCould not fit model with xgboost:\n",
                    "Formula: ",character_formula))
    }
    predicted_values <- vector(mode = "numeric",length = NROW(intervened_data))
    ivars <- all.vars(formula(character_formula))[-1]
    # remove all other variables to avoid false positive missing values
    intervened_data <- intervened_data[,ivars,with = FALSE]
    # then check for missing values
    no_missing <- !rowSums(is.na(intervened_data))
    intervened_data = intervened_data[no_missing]
    predicted_values[!no_missing] <- as.numeric(NA)
    rhs <- stats::reformulate(attr(stats::terms(stats::formula(character_formula)), "term.labels"))
    if (inherits(try(
        new_model_matrix <- Publish::specialFrame(rhs,
                                                  data = intervened_data,
                                                  specials = NULL,
                                                  na.action = na.omit)$design
    ),"try-error")){
        stop(paste0("Problem with the design matrix for the xgboost predictions.\n",
                    "A possible reason could be a variable that has no variation in the intervened_data."))
    }
    ## nd <- xgboost::xgb.DMatrix(data = new_model_matrix)
    if (inherits(try(
        predicted_values[no_missing] <- predict(fit,newdata = new_model_matrix)
    ),"try-error")) {
        stop("Xgboost prediction failed")
    }
    data.table::setattr(predicted_values,"fit",NULL)
    return(predicted_values)
}


######################################################################
### learn_xgboost.R ends here
