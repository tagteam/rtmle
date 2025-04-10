##' Method implementing the cheap subsampling method for confidence intervals
##'
##' Given a model object or a function that returns a vector of coefficients
##' and a data set,
##' this function computes confidence intervals using
##' the cheap subsampling method.
##'
##' @title Cheap subsampling bootstrap
#' @param x object of class \code{rtmle}
#' @param B Number of bootstrap samples.
#' @param M Size of the bootstrap samples. If argument \code{replace}
#'     is \code{FALSE} this defaults to \code{0.632 *
#'     NROW(x$prepared_data)} and if argument \code{replace} is
#'     \code{TRUE} this defaults to \code{NROW(x$prepared_data)}.
#' @param replace Logical. Default is \code{FALSE}. If \code{FALSE}
#'     bootstrap samples are obtained by sampling without
#'     replacement. In this case the value of argument \code{M} must
#'     be smaller than \code{NROW(x$prepared_data)}.  # If \code{TRUE}
#'     bootstrap samples are obtained by sampling with replacement. In
#'     this case argument, by default \code{M} will be set equal to
#'     \code{NROW(x$prepared_data)} # Defaults to "subsampling".
#' @param seeds A vector of random seeds for controlling the
#'     randomness while drawing the bootstrap samples.
#' @param alpha Significance level. Defaults to 0.05.
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
#' x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
#' add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
#' x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
#' protocol(x) <- list(name = "Always_A",
#'                     intervention = data.frame("A" = factor("1",levels = c("0","1"))))
#' protocol(x) <- list(name = "Never_A",
#'                     intervention = data.frame("A" = factor("0",levels = c("0","1"))))
#' prepare_data(x) <- list()
#' target(x) <- list(name = "Outcome_risk",
#'                   strategy = "additive",
#'                   estimator = "tmle",
#'                   protocols = c("Always_A","Never_A"))
#'  x <- run_rtmle(x,learner = "learn_glm",time_horizon = 1:tau)
#'  x <- cheap_bootstrap(x,B=5,time_horizon=1:tau)
#'  summary(x)
#' cb
#' }
#'
#' @export
cheap_bootstrap <- function(x,
                            B = 20,
                            M,
                            replace = FALSE,
                            seeds,
                            alpha = 0.05,
                            verbose = FALSE){
    if (missing(seeds)) seeds <- sample.int(1e+09, B)
    N <- NROW(x$prepared_data)
    for (b in 1:B){
        set.seed(seeds[[b]])
        if (replace[[1]] == TRUE) {
            inbag <- sample(1:N,size = N,replace = TRUE)
        }else{
            if(M >= N) stop(paste0("When replace = FALSE, the subsample size M (specified is M=",M,") must be lower than the sample size N (which is N=",N,")"))
            inbag <- sample(1:N,size = M,replace = FALSE)
        }
        x <- run_rtmle(x = x, verbose = verbose,
                       subsets = list(list(label = "Cheap_bootstrap",id = x$prepared_data[[x$names$id]][inbag],append = TRUE,variable = "B",value = b)))
    }
    # calculate the cheap lower and upper confidence limits
    # when there readily are bootstrap results we append to them
    if("Main"%in%names(x$estimate$Cheap_bootstrap)){
        data.table::set(x$estimate$Cheap_bootstrap,j = "Main",value = NULL)
    }
    cb <- x$estimate[["Main_analysis"]][,data.table::data.table(Target,Protocol,Time_horizon,Main = Estimate)][x$estimate$Cheap_bootstrap,on = c("Target","Protocol","Time_horizon")]
    if (replace == FALSE) {
        cheap_scale <- sqrt(M / (N - M))
    } else {
        cheap_scale <- 1
    }
    # t-distribution 
    cb[,tq := stats::qt(1 - alpha / 2, df = B)]
    cb[,cheap_variance := cumsum((Main-Estimate)^2)/seq_len(.N),by = c("Target","Protocol","Time_horizon")]
    cb[,cheap_lower := Estimate - tq * cheap_scale * cheap_variance]
    cb[,cheap_upper := Estimate + tq * cheap_scale * cheap_variance]
    x$estimate$Cheap_bootstrap <- cb[,data.table::data.table(B,
                                        Target,
                                        Protocol,
                                        Time_horizon,
                                        Target_parameter,
                                        Estimator,
                                        Estimate = Main,
                                        Bootstrap_estimate = Estimate,
                                        Bootstrap_standard_error = sqrt(cheap_variance),
                                        Bootstrap_lower = cheap_lower,
                                        Bootstrap_upper = cheap_upper)]
    x$estimate[["Main_analysis"]] <- x$estimate$Cheap_bootstrap[.N,data.table::data.table(Target,Protocol,Time_horizon,Bootstrap_lower,Bootstrap_upper)][x$estimate[["Main_analysis"]],on = c("Target","Protocol","Time_horizon")]
    # NOTE if we would return x[] instead of x then x looses its class!
    return(x)
}
