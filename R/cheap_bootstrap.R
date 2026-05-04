##' Cheap subsampling bootstrap confidence intervals
##'
##' Computes confidence intervals with the cheap subsampling bootstrap for a
##' fitted \code{rtmle} object.
##'
##' @title Cheap subsampling bootstrap
#' @param x An object of class \code{"rtmle"} with estimates produced by
#'   \code{\link{run_rtmle}}.
#' @param time_horizon The time horizon at which to calculate bootstrap
#'     confidence intervals. If this is a vector, the analysis is performed for
#'     each element.
#' @param B Number of bootstrap samples.
#' @param M Size of the bootstrap samples. The default is
#'   \code{0.632 * NROW(x$prepared_data)}.
#' @param seeds A vector of random seeds for controlling the
#'     randomness while drawing the bootstrap samples.
#' @param alpha Significance level. Defaults to 0.05.
#' @param add Logical. If \code{TRUE}, add new bootstrap results to existing
#'   bootstrap results.
#' @param verbose Logical. Passed to \code{\link{run_rtmle}} to control
#'   verbosity.
#' @param analyses Character vector specifying which analyses should receive
#'   cheap-bootstrap confidence intervals.
#' @param replace Logical. Whether to sample with replacement. The default is
#'   \code{FALSE} and should usually not be changed.
#' @param ... Additional arguments passed to \code{run_rtmle}.
##' @return The modified \code{rtmle} object.
##' @seealso \code{\link{run_rtmle}}, \code{\link{summary.rtmle}}
##' @author Johan Sebastian Ohlendorff
##'     \email{johan.ohlendorff@@sund.ku.dk} and Thomas A Gerds
##'     \email{tag@@biostat.ku.dk}
#' @examples
#' \dontrun{
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
                            analyses = NULL,
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
    if (is.null(analyses)){
        analyses <- names(x$estimate)
    } else {
        ## Make sure that the specified subsets exist
        analyses <- intersect(analyses, names(x$estimate))
        if (length(analyses) == 0){
            stop("None of the specified subsets in 'analyses' exist in the object.")
        }
    }
    analyses <- setdiff(analyses,"Cheap_bootstrap")
    Target_parameter <- Main_estimate <- Bootstrap_estimate <- Bootstrap_lower <- Bootstrap_upper <- Main <- cheap_lower <- cheap_upper <- cheap_variance <- Estimate <- Estimator <- tq <- Time_horizon <- Protocol <- Target <- NULL
    for (v in analyses){
        if (verbose) message(paste0("Cheap bootstrap for '",v,"'"))
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
            x <- run_rtmle(x = x,refit = TRUE,verbose = verbose,time_horizon = time_horizon,subsets = list(list(label = "CB",id = x$prepared_data[[x$names$id]][inbag],append = TRUE,variable = "B",level = b)),learner = x$unparsed_learner,...)
        }
        x$estimate$Cheap_bootstrap[[v]] <- x$estimate$CB
        x$estimate$CB <- NULL
      # calculate the cheap lower and upper confidence limits
      # append to existing bootstrap results when present
      if("Main"%in%names(x$estimate$Cheap_bootstrap[[v]])){
        data.table::set(x$estimate$Cheap_bootstrap[[v]],j = "Main",value = NULL)
      }
      data.table::setnames(x$estimate$Cheap_bootstrap[[v]],"Estimate","Bootstrap_estimate")
      cb <- x$estimate[[v]][,list(Target,Protocol,Time_horizon,Main_estimate = Estimate)][x$estimate$Cheap_bootstrap[[v]],on = c("Target","Protocol","Time_horizon")]
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
      x$estimate$Cheap_bootstrap[[v]] <- cb[,list(B,
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
      cb <- cb[B == max_B,list(Bootstrap_standard_error = sqrt(cheap_variance),
                                                 Bootstrap_lower = pmax(0,cheap_lower),
                                                 Bootstrap_upper = pmin(1,cheap_upper)),by = c("Target","Protocol","Time_horizon")]
      x$estimate[[v]] <- cb[x$estimate[[v]],on = c("Target","Protocol","Time_horizon")]
    }
    # NOTE if we would return x[] instead of x then x looses its class!
    return(x)
}
