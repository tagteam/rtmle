### learn_xgboost.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 28 2024 (09:26) 
## Version: 
## Last-Updated: nov 24 2025 (14:28) 
##           By: Thomas Alexander Gerds
##     Update #: 105
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Learning nuisance parameter models for TMLE and predicting
##' probabilities in intervened data based on \code{\link[xgboost]{xgboost}}
##'
##' Hyperparameters are set via ...
##' @title Nuisance parameter learner based on \code{\link[xgboost]{xgboost}}
##' @param character_formula Formula for nuisance parameter as a character
##' @param data Data for learning 
##' @param intervened_data Data for prediction 
##' @param ... Additional arguments for the learning phase passed to \code{\link[xgboost]{xgboost}}. These can include hyperparameters.
##' @return A vector of predicted probabilities which has the fit as an attribute.  
##' @seealso \code{link{superlearn}}, \code{link{learn_glm}}, \code{link{learn_glmnet}}
##' @examples
#' set.seed(17)
#' tau <- 3
#' ld <- simulate_long_data(n = 391,number_visits = 20,
#'                          beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),
#'                          register_format = TRUE)
#' x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",
#'                 name_censoring = "Censored",censored_label = "censored")
#' x <- add_long_data(x,
#'                    outcome_data=ld$outcome_data,
#'                    censored_data=ld$censored_data,
#'                    competing_data=ld$competing_data,
#'                    timevar_data=ld$timevar_data)
#' x <- add_baseline_data(x,data=ld$baseline_data)
#' x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
#' x <- protocol(x,name = "Always_A",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1"))))
#' x <- protocol(x,name = "Never_A",
#'                     intervention = data.frame("A" = factor("0",levels = c("0","1"))))
#' x <- prepare_data(x)
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
                                                                  eta = 1.0,
                                                                  nrounds = 2L,
                                                                  subsample = 1.0,
                                                                  lambda = 1,
                                                                  verbose = 0,
                                                                  objective = "reg:squarederror")),
                                     forced=NULL,
                                     verbose=TRUE)
    # extract the data
    # FIXME: should not let NA's pass until here
    sf <- Publish::specialFrame(stats::formula(character_formula),
                                data = data,
                                specials = NULL,
                                na.action = na.omit)
    Y <- sf$response[[1]]
    if (is.factor(Y)){
        # FIXME: the reference level should be controlled better
        Y <- as.numeric(Y == levels(Y)[[2]])
    }
    d <- xgboost::xgb.DMatrix(sf$design, label = Y)
    if (!inherits(try(
             fit <- do.call(xgboost::xgboost,c(list(data = d),control$xgboost))
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
    no_missing <- !(apply(intervened_data,1,function(x)any(is.na(x))))
    predicted_values[!no_missing] <- as.numeric(NA)
    rhs <- stats::reformulate(attr(stats::terms(stats::formula(character_formula)), "term.labels"))
    new_model_matrix <- Publish::specialFrame(rhs,
                                              data = intervened_data[no_missing],
                                              specials = NULL,
                                              na.action = na.omit)$design
    nd <- xgboost::xgb.DMatrix(new_model_matrix)
    if (inherits(try(
        predicted_values[no_missing] <- predict(fit,newdata = nd)
    ),"try-error")) {
        stop("Xgboost prediction failed")
    }
    data.table::setattr(predicted_values,"fit",NULL)
    return(predicted_values)
}


######################################################################
### learn_xgboost.R ends here
