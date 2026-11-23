##' Cheap subsampling bootstrap confidence intervals
##'
##' Computes confidence intervals with the cheap subsampling bootstrap for a
##' fitted \code{rtmle} object. All estimates selected for bootstrapping must use the same estimation method.
##'
##' @title Cheap subsampling bootstrap
#' @param x An object of class \code{"rtmle"} with estimates produced by
#'   \code{\link{run_rtmle}}.
#' @param time_horizon The time horizon at which to calculate bootstrap
#'     confidence intervals. If this is a vector, the analysis is performed for
#'     each element. Requested horizons must have fitted estimates.
#' @param B Number of bootstrap samples.
#' @param M Size of the bootstrap samples. The default is
#'   \code{0.632 * NROW(x$prepared_data)}.
#' @param seeds A vector of random seeds for controlling the
#'     randomness while drawing the bootstrap samples.
#' @param alpha Significance level. Defaults to 0.05.
#' @param add Logical. If \code{TRUE}, add new bootstrap results to existing
#'   bootstrap results for the same estimates and resampling settings.
#' @param verbose Logical. Passed to \code{\link{run_rtmle}} to control
#'   verbosity.
#' @param replace Logical. Whether to sample with replacement. The default is
#'   \code{FALSE} and should usually not be changed.
#' @param targets Character vector specifying which fitted targets should receive
#'   cheap-bootstrap confidence intervals. Defaults to all fitted targets.
#' @param ... Additional arguments passed to \code{run_rtmle}.
##' @return The modified \code{rtmle} object.
##' @seealso \code{\link{run_rtmle}}, \code{\link{summary.rtmle}}
##' @author Johan Sebastian Ohlendorff
##'     \email{johan.ohlendorff@@sund.ku.dk} and Thomas A Gerds
##'     \email{tag@@biostat.ku.dk}
#' @examples
#' \dontrun{
#' data(rtmle_object)
#' rtmle_object <- cheap_bootstrap(rtmle_object,
#'                                  B=5,
#'                                  time_horizon=1,
#'                                  seeds=17)
#' summary(rtmle_object)
#' }
#'
#' @export
cheap_bootstrap <- function(x, ...) {
    UseMethod("cheap_bootstrap")
}

#' @rdname cheap_bootstrap
#' @export
cheap_bootstrap.rtmle <- function(x,
                                  time_horizon,
                                  B = 20,
                                  M,
                                  seeds,
                                  alpha = 0.05,
                                  add = TRUE,
                                  verbose = FALSE,
                                  replace = FALSE,
                                  targets = NULL,
                                  ...) {
    Target <- Regime <- Time_horizon <- Estimator <- Target_parameter <- NULL
    Estimate <- Main_Estimate <- Bootstrap_estimate <- Bootstrap_standard_error <- NULL

    main <- x$estimate$Main_analysis
    if (NROW(main) == 0L) {
        stop("Run run_rtmle() before calling cheap_bootstrap().")
    }
    if (replace) {
        message("Argument 'replace' should not be changed from its default value FALSE for cheap bootstrap. \n Only valid when not using a super learner.")
    }
    fit_args <- list(...)
    if ("analyses" %in% names(fit_args)) {
        stop("Call cheap_bootstrap() on the fitted object for the analysis of interest.")
    }
    fitted <- main[is.null(targets) | as.character(Target) %in% targets,
                   list(Target, Regime, Time_horizon, Estimator, Main_Estimate = Estimate)]
    if (missing(time_horizon)) {
        time_horizon <- sort(unique(fitted$Time_horizon))
    } else if (NROW(fitted) > 0L) {
        missing_horizons <- setdiff(time_horizon, fitted$Time_horizon)
        if (length(missing_horizons) > 0L) {
            stop("No fitted estimates for time horizon(s): ",
                 paste(missing_horizons, collapse = ", "))
        }
        fitted <- fitted[Time_horizon %in% time_horizon]
    }
    fitted_targets <- unique(as.character(fitted$Target))
    missing_targets <- setdiff(targets, fitted_targets)
    if (length(missing_targets) > 0L) {
        stop("No estimates at the selected time horizons for target(s): ",
             paste(missing_targets, collapse = ", "))
    }
    if (is.null(targets)) targets <- fitted_targets
    estimator <- unique(as.character(fitted$Estimator))
    if (length(estimator) != 1L) {
        stop("Select fitted targets and horizons for one estimator.")
    }

    N <- NROW(x$prepared_data)
    M <- if (missing(M)) round(0.632 * N) else round(M)
    stopifnot(length(B) == 1L, B > 0, B == floor(B), length(M) == 1L, M > 0)
    if (!replace && M >= N) {
        stop("The subsample size M must be smaller than the sample size N.")
    }
    cheap_scale <- if (replace) sqrt(M / N) else sqrt(M / (N - M))
    resampling <- list(
        N = N, M = M, replace = replace,
        fit_args = fit_args[intersect(c("seed", names(x$tuning_parameters)), names(fit_args))]
    )

    keys <- c("Target", "Regime", "Time_horizon", "Estimator")
    draw_columns <- c("B", keys, "Target_parameter", "Bootstrap_estimate")
    previous <- if (add) x$estimate$Cheap_bootstrap$Main_analysis
    offset <- 0L
    if (NROW(previous) > 0L) {
        previous_resampling <- attr(previous, "resampling", exact = TRUE)
        if ((!is.null(previous_resampling) &&
             !identical(previous_resampling, resampling[names(previous_resampling)])) ||
            !isTRUE(all.equal(attr(previous, "cheap_scale"), cheap_scale)) ||
            !data.table::fsetequal(unique(previous[, names(fitted), with = FALSE]), fitted)) {
            stop("Use add = FALSE to change the fitted estimates or the resampling settings.")
        }
        offset <- max(previous$B)
        previous <- previous[, draw_columns, with = FALSE]
    }
    total <- offset + B
    if (missing(seeds)) {
        seeds <- sample.int(1e9, total)
    } else if (length(seeds) == 1L) {
        set.seed(seeds)
        seeds <- sample.int(1e9, total)
    }
    stopifnot(length(seeds) == total)

    draws <- lapply(offset + seq_len(B), function(b) {
        set.seed(seeds[[b]])
        rows <- sample.int(N, M, replace = replace)
        xb <- subset_rtmle_rows(x, rows)
        xb$prepared_data[[x$names$id]] <- seq_along(rows)
        if (NROW(xb$followup) > 0L) xb$followup[[x$names$id]] <- seq_along(rows)
        xb <- run_rtmle(xb, targets = targets,
                        learner = x$unparsed_learner, estimator = estimator,
                        time_horizon = time_horizon,
                        verbose = verbose, ...)
        xb$estimate$Main_analysis[, list(
            B = b, Target, Regime, Time_horizon, Estimator, Target_parameter,
            Bootstrap_estimate = Estimate
        )]
    })

    cb <- fitted[data.table::rbindlist(c(list(previous), draws)), on = keys, nomatch = 0L]
    tq <- stats::qt(1 - alpha / 2, df = total)
    cb[, Bootstrap_standard_error :=
           cheap_scale * sqrt(cumsum((Main_Estimate - Bootstrap_estimate)^2) / seq_len(.N)),
       by = keys]
    cb[, `:=`(
        Bootstrap_lower = pmax(0, Main_Estimate - tq * Bootstrap_standard_error),
        Bootstrap_upper = pmin(1, Main_Estimate + tq * Bootstrap_standard_error)
    )]
    data.table::setcolorder(cb, c("B", keys, "Target_parameter", "Main_Estimate",
                                "Bootstrap_estimate"))
    data.table::setattr(cb, "cheap_scale", cheap_scale)
    data.table::setattr(cb, "tq", tq)
    data.table::setattr(cb, "resampling", resampling)
    x$estimate$Cheap_bootstrap <- list(Main_analysis = cb)

    intervals <- c("Bootstrap_standard_error", "Bootstrap_lower", "Bootstrap_upper")
    main <- main[, setdiff(names(main), intervals), with = FALSE]
    last_draw <- cb[cb$B == total, c(keys, intervals), with = FALSE]
    x$estimate$Main_analysis <- last_draw[main, on = keys]
    x
}
