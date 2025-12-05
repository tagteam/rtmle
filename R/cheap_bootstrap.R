##' Method implementing the cheap subsampling method for confidence intervals
##'
##' Given a model object or a function that returns a vector of coefficients
##' and a data set,
##' this function computes confidence intervals using
##' the cheap subsampling method.
##'
##' @title Cheap subsampling bootstrap
#' @param x object of class \code{rtmle}
#' @param time_horizon  The time horizon at which to calculate
#'     the bootstrap confidence intervals. If it is a vector the analysis will be performed for
#'     each element of the vector.
#' @param B Number of bootstrap samples.
#' @param M Size of the bootstrap samples. The defaults is \code{0.632 * NROW(x$prepared_data)}.
#' @param seeds A vector of random seeds for controlling the
#'     randomness while drawing the bootstrap samples.
#' @param alpha Significance level. Defaults to 0.05.
#' @param add Logical. If \code{TRUE} add new bootstrap results to the existing bootstrap results
#' @param verbose Logical. Passed to run_rtmle to control verbosity.
##' @return The modified object
##' @author Johan Sebastian Ohlendorff
##'     \email{johan.ohlendorff@@sund.ku.dk} and Thomas A Gerds
##'     \email{tag@@biostat.ku.dk}
#' @examples
#' \dontrun{
#' set.seed(17)
#' tau <- 3
#' ld <- simulate_long_data(n = 100,number_visits = 20,
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
#'  x <- model_formula(x)
#'  x <- run_rtmle(x,learner = "learn_glm",time_horizon = 1:tau)
#'  x <- cheap_bootstrap(x,B=5,time_horizon=1:tau)
#'  summary(x)
#' cb
#' }
#'
#' @export
cheap_bootstrap <- function(x,
                            time_horizon,
                            B = 20,
                            M,
                            seeds,
                            alpha = 0.05,
                            add = TRUE,
                            verbose = FALSE){
    if (length(x$estimate$Main_analysis) == 0){
        stop("The object lacks the main analysis. Please apply run_rtmle() before calling cheap_bootstrap.")
    }
    Target_parameter <- Main_estimate <- Bootstrap_estimate <- Bootstrap_lower <- Bootstrap_upper <- Main <- cheap_lower <- cheap_upper <- cheap_variance <- Estimate <- Estimator <- tq <- Time_horizon <- Protocol <- Target <- NULL
    # add to existing bootstrap results if any
    if (add[[1]] == TRUE && length(x$estimate$Cheap_bootstrap)>0){
        B_offset <- max(x$estimate$Cheap_bootstrap$B)
        x$estimate$Cheap_bootstrap[,"P_value" := NA]
        x$estimate$Cheap_bootstrap[, "Main_Estimate":= NULL]
        setnames(
          x$estimate$Cheap_bootstrap,
          c(
            "Bootstrap_estimate",
            "Bootstrap_standard_error",
            "Bootstrap_lower",
            "Bootstrap_upper"
          ),
          c("Estimate", 
            "Standard_error",
            "Lower",
            "Upper"
          )
        )
    }else{
        x$estimate$Cheap_bootstrap <- NULL
        B_offset <- 0
    }
    if (missing(seeds)) seeds <- sample.int(1e+09, B+B_offset)
    N <- NROW(x$prepared_data)
    for (b in (1:B)+B_offset){
        set.seed(seeds[[b]])
        if(M >= N) stop(paste0("The subsample size M (specified is M=",M,") must be lower than the sample size N (which is N=",N,")"))
        inbag <- sample(1:N,size = M,replace = FALSE)
        if (missing(time_horizon)){
            time_horizon <- sort(unique(x$estimate$Main_analysis$Time_horizon))
        }
        x <- run_rtmle(x = x,
                       verbose = verbose,
                       time_horizon = time_horizon,
                       subsets = list(list(label = "Cheap_bootstrap",
                                           id = x$prepared_data[[x$names$id]][inbag],
                                           append = TRUE,
                                           variable = "B",
                                           level = b)))
    }
    # calculate the cheap lower and upper confidence limits
    # when there readily are bootstrap results we append to them
    if("Main"%in%names(x$estimate$Cheap_bootstrap)){
        data.table::set(x$estimate$Cheap_bootstrap,j = "Main",value = NULL)
    }
    data.table::setnames(x$estimate$Cheap_bootstrap,"Estimate","Bootstrap_estimate")
    cb <- x$estimate[["Main_analysis"]][,data.table::data.table(Target,Protocol,Time_horizon,Main_estimate = Estimate)][x$estimate$Cheap_bootstrap,on = c("Target","Protocol","Time_horizon")]
    cheap_scale <- sqrt(M / (N - M))
    # t-distribution 
    cb[,tq := stats::qt(1 - alpha / 2, df = B)]
    cb[,cheap_variance := cumsum((Main_estimate-Bootstrap_estimate)^2)/seq_len(.N),by = c("Target","Protocol","Time_horizon")]
    cb[,cheap_lower := Main_estimate - tq * cheap_scale * sqrt(cheap_variance)]
    cb[,cheap_upper := Main_estimate + tq * cheap_scale * sqrt(cheap_variance)]
    x$estimate$Cheap_bootstrap <- cb[,data.table::data.table(B,
                                                             Target,
                                                             Protocol,
                                                             Time_horizon,
                                                             Target_parameter,
                                                             Estimator,
                                                             Main_Estimate = Main_estimate,
                                                             Bootstrap_estimate = Bootstrap_estimate,
                                                             Bootstrap_standard_error = sqrt(cheap_variance),
                                                             Bootstrap_lower = pmax(0,cheap_lower),
                                                             Bootstrap_upper = pmin(1,cheap_upper))]
    max_B <- B+B_offset
    data.table::setattr(x$estimate$Cheap_bootstrap,
                        "cheap_scale",cheap_scale)
    data.table::setattr(x$estimate$Cheap_bootstrap,
                        "tq",unique(cb[B == max_B,tq]))
    # merge the bootstrap confidence limits
    # first remove what is there already
    if ("Bootstrap_lower" %in% names(x$estimate[["Main_analysis"]])){
        for (v in c("Bootstrap_lower","Bootstrap_upper","Bootstrap_standard_error")){
            data.table::set(x$estimate[["Main_analysis"]],j = v,value = NULL)
        }
    }
    cb <- cb[B == max_B,data.table::data.table(Bootstrap_standard_error = sqrt(cheap_variance),
                                               Bootstrap_lower = pmax(0,cheap_lower),
                                               Bootstrap_upper = pmin(1,cheap_upper)),by = c("Target","Protocol","Time_horizon")]
    x$estimate[["Main_analysis"]] <- cb[x$estimate[["Main_analysis"]],on = c("Target","Protocol","Time_horizon")]
    # NOTE if we would return x[] instead of x then x looses its class!
    return(x)
}
