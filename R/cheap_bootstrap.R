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
#' @param which_subsets Character vector specifying which subsets to perform the cheap bootstrap on.
#' @param replace Logical. Whether to sample with replacement. Default is FALSE. Should not be changed.
#' @param ... Additional arguments passed to \code{run_rtmle}.
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
#'  x <- cheap_bootstrap(x,B=5,time_horizon=1:tau, M = 0.632*NROW(x$prepared_data))
#'  summary(x)
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
                            verbose = FALSE,
                            which_subsets = NULL,
                            replace = FALSE,
                            ...){
    if (is.null(x$estimate)){
        stop("The object lacks an analysis. Please apply run_rtmle() before calling cheap_bootstrap.")
    }
    if (!("Cheap_bootstrap") %in% names(x$estimate)){
      x$estimate$Cheap_bootstrap <- list()
    }
    if (replace){
      message("Argument 'replace' should not be changed from its default value FALSE for cheap bootstrap. \n Only valid when not using a super learner.")
    } 
    if (is.null(which_subsets)){
      which_subsets <- names(x$estimate)
    } else {
      ## Make sure that the specified subsets exist
      which_subsets <- intersect(which_subsets, names(x$estimate))
      if (length(which_subsets) == 0){
        stop("None of the specified subsets in 'which_subsets' exist in the object.")
      }
    }
    which_subsets <- setdiff(which_subsets,"Cheap_bootstrap")
    
    Target_parameter <- Main_estimate <- Bootstrap_estimate <- Bootstrap_lower <- Bootstrap_upper <- Main <- cheap_lower <- cheap_upper <- cheap_variance <- Estimate <- Estimator <- tq <- Time_horizon <- Protocol <- Target <- NULL
    for (v in which_subsets){
      message(paste0("Cheap bootstrap for '",v,"'"))
      # add to existing bootstrap results if any
      if (add[[1]] == TRUE && length(x$estimate$Cheap_bootstrap[[v]])>0){
        B_offset <- max(x$estimate$Cheap_bootstrap[[v]]$B)
        x$estimate$Cheap_bootstrap[[v]][,"P_value" := NA]
        x$estimate$Cheap_bootstrap[[v]][, "Main_Estimate":= NULL]
        setnames(
          x$estimate$Cheap_bootstrap[[v]],
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
        x$estimate$Cheap_bootstrap[[v]] <- NULL
        B_offset <- 0
      }
      if (missing(seeds)) seeds <- sample.int(1e+09, B+B_offset)
      N <- NROW(x$prepared_data)
      for (b in (1:B)+B_offset){
        set.seed(seeds[[b]])
        if(M >= N && !replace) stop(paste0("The subsample size M (specified is M=",M,") must be lower than the sample size N (which is N=",N,")"))
        inbag <- sample(1:N,size = M,replace = replace)
        if (missing(time_horizon)){
          time_horizon <- sort(unique(x$estimate[[v]]$Time_horizon))
        }
        x <- run_rtmle(x = x,
                       verbose = verbose,
                       time_horizon = time_horizon,
                       subsets = list(list(label = "CB",
                                           id = x$prepared_data[[x$names$id]][inbag],
                                           append = TRUE,
                                           variable = "B",
                                           level = b)),
                       learner = x$learner,
                       ...)
      }
      x$estimate$Cheap_bootstrap[[v]] <- x$estimate$CB
      x$estimate$CB <- NULL

      # calculate the cheap lower and upper confidence limits
      # when there readily are bootstrap results we append to them
      if("Main"%in%names(x$estimate$Cheap_bootstrap[[v]])){
        data.table::set(x$estimate$Cheap_bootstrap[[v]],j = "Main",value = NULL)
      }
      data.table::setnames(x$estimate$Cheap_bootstrap[[v]],"Estimate","Bootstrap_estimate")
      cb <- x$estimate[[v]][,data.table::data.table(Target,Protocol,Time_horizon,Main_estimate = Estimate)][x$estimate$Cheap_bootstrap[[v]],on = c("Target","Protocol","Time_horizon")]
      if (replace){
        cheap_scale <- sqrt(M / N)
      } else {
        cheap_scale <- sqrt(M / (N - M))
      }
      # t-distribution 
      cb[,tq := stats::qt(1 - alpha / 2, df = B)]
      cb[,cheap_variance := cumsum((Main_estimate-Bootstrap_estimate)^2)/seq_len(.N),by = c("Target","Protocol","Time_horizon")]
      cb[,cheap_lower := Main_estimate - tq * cheap_scale * sqrt(cheap_variance)]
      cb[,cheap_upper := Main_estimate + tq * cheap_scale * sqrt(cheap_variance)]
      x$estimate$Cheap_bootstrap[[v]] <- cb[,data.table::data.table(B,
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
      data.table::setattr(x$estimate$Cheap_bootstrap[[v]],
                          "cheap_scale",cheap_scale)
      data.table::setattr(x$estimate$Cheap_bootstrap[[v]],
                          "tq",unique(cb[B == max_B,tq]))
      # merge the bootstrap confidence limits
      # first remove what is there already
      if ("Bootstrap_lower" %in% names(x$estimate[[v]])){
        for (u in c("Bootstrap_lower","Bootstrap_upper","Bootstrap_standard_error")){
          data.table::set(x$estimate[[v]],j = u,value = NULL)
        }
      }
      cb <- cb[B == max_B,data.table::data.table(Bootstrap_standard_error = sqrt(cheap_variance),
                                                 Bootstrap_lower = pmax(0,cheap_lower),
                                                 Bootstrap_upper = pmin(1,cheap_upper)),by = c("Target","Protocol","Time_horizon")]
      x$estimate[[v]] <- cb[x$estimate[[v]],on = c("Target","Protocol","Time_horizon")]
    }
    # NOTE if we would return x[] instead of x then x looses its class!
    return(x)
}
