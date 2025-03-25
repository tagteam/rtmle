## retrieve cheap subsampling confidence interval
cheap_subsampling_ci <- function(est, boot_est, size, n_val, alpha) {
  b_val <- length(boot_est)
  s_val <- sqrt(mean((est - boot_est)^2))
  tq <- stats::qt(1 - alpha / 2, df = b_val)
  data.frame(
    estimate = est,
    cheap_lower = est - tq * sqrt((size) / (n_val - size)) * s_val,
    cheap_upper = est + tq * sqrt((size) / (n_val - size)) * s_val
  )
}

## get data and retrieve function from call object in fun.
get_data_and_call <- function(fun, data = NULL, parent_env) {
  if (!inherits(fun, "function")) {
    ## tryCatch to retrieve fun$call
    message("fun is not a function, trying to retrieve call object from fun")
    call <- fun$call
    if (is.null(call)) {
      stop("fun does not appear to have a call")
    }
    # get data from call in parent environment
    if (is.null(data)) {
      tryCatch(
        {
          data <- eval(call$data, parent_env)
        },
        error = function(e) {
          stop("data not found in environment or missing from call object")
        }
      )
    }
    ## make call object into string
    fun <- function(d) {
      call$data <- force(d)
      do.call("coef", list(eval(call, envir = parent_env)), envir = parent_env)
    }
  }

  if (!inherits(data, "data.frame")) {
    stop("Data needs to be a data frame")
  }
  n_val <- nrow(data)
  list(
    data = data,
    n_val = n_val,
    fun = fun
  )
}

##' Method implementing the cheap subsampling method for confidence intervals
##'
##' Given a model object or a function that returns a vector of coefficients
##' and a data set,
##' this function computes confidence intervals using
##' the cheap subsampling method.
##'
##' @title Cheap subsampling
##' @param fun A function that returns a vector of coefficients or a
##'     model object which saves the call. In the second case, a coef
##'     needs to be defined for 'class(fun)'.
##' @param data Data set to be used for the computation, if applicable.
##' @param b Number of bootstrap samples.
##' @param size Subsample size. Defaults to 0.632 * nrow(data).
##' @param alpha Significance level. Defaults to 0.05.
##' @param type Character. Type of bootstrap method.
##' Can be either "subsampling" or "non_parametric".
##' Defaults to "subsampling".
##' @param parallel_args List of arguments for parallel computation.
##' The first element parallelize should be a logical.
##' Should the bootstrap samples be computed in parallel?
##' Defaults to FALSE.
##' The second element cores should be an integer.
##' Number of cores to use for parallel computation,
##' if parallelize = TRUE. Defaults to 2.
##' @param adapt_args List of arguments for adaptive computation.
##' The first element adapt should be a logical.
##' Should the number of bootstrap samples be adaptively chosen?
##' Defaults to FALSE.
##' The second element precision should be a numeric.
##' Desired precision of the confidence interval(s). Defaults to 0.1.
##' The third element max_b should be an integer.
##' Maximum number of bootstrap samples to be taken. Defaults to 200.
##' @param additional_args List of additional arguments.
##' The first element keep_boot_estimates should be a logical.
##' Should the bootstrap estimates be kept? Defaults to FALSE.
##' The second element verbose should be a logical.
##' Should a progress bar be displayed or if adapt = TRUE,
##' the current precision printed.
##' @return An object of class "cheap_bootstrap" containing
##' the point estimates and confidence intervals.
##' @export
##' @examples
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' fun <- function(data) {
##'   coef(lm(formula = Postwt ~ Prewt + Treat + offset(Prewt), data = data))
##' }
##' cs <- cheap_bootstrap(fun = fun, b = 20, data = anorexia)
##' cs
##' summary(cs)
##'
##' ## example with a model object
##' set.seed(123)
##' x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
##' cs2 <- cheap_bootstrap(fun = x, b = 20)
##' summary(cs2)
##'
##' ## example with adaptive method
##' set.seed(123)
##' fun <- function(data) {
##'   coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = data))
##' }
##' cs <- cheap_bootstrap(fun = fun,
##'                       data = anorexia,
##'                       adapt_args = list(adapt = TRUE,
##'                                         precision = 0.1,
##'                                         max_b = 200))
##' cs
##' summary(cs)
##'
##' ## example with a function call and parallel computation
##' \dontrun{
##' set.seed(123)
##' ## note the function needs to load the packages needed
##' ## for the computation if parallelized
##' x <- function(data) {
##'   coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = data))
##' }
##' cs3 <- cheap_bootstrap(fun = x,
##'                        b = 20,
##'                        data = anorexia,
##'                        parallel_args = list(parallelize = TRUE,
##'                                             cores = 2))
##' summary(cs3)
##'
##' ## example with custom function returning coefficients
##' ## example from ATE-package
##' library(survival)
##' library(riskRegression)
##' set.seed(10)
##' #### Survival settings  ####
##' #### ATE with Cox model ####
##'
##' ## generate data
##'
##' n <- 100
##' dt <- sampleData(n, outcome="survival")
##' dt$time <- round(dt$time, 1)
##' dt$X1 <- factor(rbinom(n,
##'                        prob = c(0.3,0.4),
##'                        size = 2),
##'                  labels = paste0("T",0:2))
##'
##' ## fitter function which returns a named vector
##' ate_fit_fun <- function(d) {
##'   ## estimate the Cox model
##'   fit <- coxph(formula = Surv(time,event)~ X1+X2,data=d,y=TRUE,x=TRUE)
##'
##'   ## main fitter function; ignore output
##'   invisible(capture.output(ate_fit <- summary(ate(fit,
##'                                                   data = d,
##'                                                   treatment = "X1",
##'                                                   times = 5:8,
##'                                                   se = FALSE))))
##'   ## extract the point estimates for risk difference
##'   res <- ate_fit$diffRisk$estimate
##'   ## name the point estimates
##'   names(res) <- paste0("ATE ",
##'                        ate_fit$diffRisk$A,
##'                        "-",
##'                        ate_fit$diffRisk$B,
##'                        " (t = ",
##'                        ate_fit$diffRisk$time,")")
##'   res
##' }
##' set.seed(102)
##' cs4 <- cheap_bootstrap(ate_fit_fun,
##'                        b = 5,
##'                        data = dt)
##' summary(cs4)
##'
##' ## example from riskRegression with no existing coef function
##' set.seed(18)
##' learndat <- sampleData(200,outcome="binary")
##' testdat <- sampleData(1200,outcome="binary")
##'
##' ## score logistic regression models
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
##' lr2 = glm(Y~X3+X5,data=learndat,family=binomial)
##' z<-Score(list("LR(X1+X2+X7+X9)"=lr1,
##'          "LR(X3+X5)"=lr2),
##'          formula=Y~1,
##'          data=testdat,
##'          metrics = "AUC")
##' coef.Score <- function(x) {
##'   res <- x$AUC$score$AUC
##'   names(res) <- paste0("AUC (model = ", x$AUC$score$model, ")")
##'   res
##' }
##' set.seed(102)
##' cs5 <- cheap_bootstrap(z, b = 20)
##' summary(cs5)
##' }
cheap_bootstrap <- function(fun,
                            data = NULL,
                            b = 20,
                            size = NULL,
                            alpha = 0.05,
                            type = "subsampling",
                            parallel_args = list(
                              parallelize = FALSE,
                              cores = 2
                            ),
                            adapt_args = list(
                              adapt = FALSE,
                              precision = 0.1,
                              max_b = 200
                            ),
                            additional_args = list(
                              keep_boot_estimates = TRUE,
                              verbose = TRUE
                            )) {
  arg_check <- c(
    "b",
    "size",
    "alpha",
    "type",
    "parallelize",
    "cores",
    "adapt",
    "precision",
    "max_b",
    "keep_boot_estimates",
    "verbose"
  )
  ## declare variables from lists
  parallelize <- parallel_args$parallelize
  cores <- parallel_args$cores
  adapt <- adapt_args$adapt
  precision <- adapt_args$precision
  max_b <- adapt_args$max_b
  keep_boot_estimates <- additional_args$keep_boot_estimates
  verbose <- additional_args$verbose


  data_and_call <- get_data_and_call(fun, data = data, parent.frame())
  fun <- data_and_call$fun
  data <- data_and_call$data
  n_val <- data_and_call$n_val

  if (is.null(size) && type == "subsampling") {
    size <- round(0.632 * n_val)
  }
  if (type == "non_parametric") {
    size <- n_val
    ## cheap_bootstrap formula defaults
    ## to m = n/2 for non-parametric bootstrap
    size_cb <- 1 / 2 * n_val
  } else {
    size_cb <- size
  }
  est <- tryCatch(
    {
      fun(data)
    },
    error = function(e) {
      stop(
        "function fun/derived from fun failed with error: ",
        conditionMessage(e),
        ".\n Did you specify fun as a function or as a model object with a corresponding call and coef function?"
      )
    }
  )

  ## check that est is a named vector
  if (!is.vector(est) || !is.numeric(est) || is.null(names(est))) {
    stop("function fun/derived from fun did not return a named vector of coefficients")
  }

  ## check that all of the arguments are of length 1; except for data and fun
  for (var_name in arg_check) {
    if (length(get(eval(var_name))) != 1) {
      stop(paste(var_name, "must be of length 1"))
    }
  }

  if (adapt) {
    b <- 1

    ## initialize empty data frame to store bootstrap replications
    boot_est <- list()
    res <- list()
    width_prev <- Inf
    while (b <= max_b) {
      ## compute bootstrap estimates
      tryCatch(
        {
          curr_boot <- fun(data[
            sample(1:n_val,
              size,
              replace = (type == "non_parametric")
            ), ,
            drop = FALSE
          ])
        },
        error = function(e) {
          stop(
            "Bootstrap computation failed with error: ",
            conditionMessage(e),
            " for bootstrap iteration b = ",
            b
          )
        }
      )
      ## apply cheap_subsampling_ci for each row of boot_est, est
      tryCatch(
        {
          boot_est <- rbind(as.data.frame(t(curr_boot)), boot_est)
          res[[b]] <- do.call(
            "rbind",
            lapply(seq_len(length(est)), function(i) {
              cheap_subsampling_ci(est[i], boot_est[, i], size_cb, n_val, alpha)
            })
          )
        },
        error = function(e) {
          stop(
            "Computation of confidence intervals failed with error: ",
            conditionMessage(e),
            ". Does your function return a vector of coefficients? ",
            "NOTE: for bootstrap iteration b = ",
            b
          )
        }
      )
      ## calculate width for current iteration
      width_curr <- res[[b]]$cheap_upper - res[[b]]$cheap_lower
      ## print current precision
      if (verbose) {
        print(abs(width_curr - width_prev))
      }
      if (all(abs(width_curr - width_prev) < precision)) {
        break
      } else {
        b <- b + 1
        width_prev <- width_curr
      }
    }
    res <- res[[length(res)]]
    if (b > max_b) {
      b <- max_b
    }
  } else if (parallelize) {
    requireNamespace("parallel")
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c("fun", "b", "size", "alpha", "data", "type"),
      envir = environment()
    )
    tryCatch(
      {
        boot_est <- do.call(
          "rbind",
          parallel::parLapply(cl, seq_len(b), function(i) {
            fun(data[
              sample(1:n_val,
                size,
                replace = (type == "non_parametric")
              ), ,
              drop = FALSE
            ])
          })
        )
        parallel::stopCluster(cl)
      },
      error = function(e) {
        parallel::stopCluster(cl)
        stop(
          "Bootstrap computation failed with error: ",
          conditionMessage(e)
        )
      }
    )
  } else {
    if (verbose) {
      pb <- utils::txtProgressBar(
        min = 1,
        max = b,
        style = 3,
        width = 50,
        char = "="
      )
    }
    tryCatch(
      {
        boot_est <- do.call("rbind", lapply(seq_len(b), function(i) {
          if (verbose) {
            utils::setTxtProgressBar(pb, i)
          }
          fun(data[
            sample(1:n_val,
              size,
              replace = (type == "non_parametric")
            ), ,
            drop = FALSE
          ])
        }))
      },
      error = function(e) {
        stop(
          "Bootstrap computation failed with error: ",
          conditionMessage(e)
        )
      }
    )
  }
  if (verbose) {
    cat("\n")
  }

  if (!adapt) {
    tryCatch(
      {
        ## apply cheap_subsampling_ci for each row of boot_est, est
        res <- do.call("rbind", lapply(seq_len(length(est)), function(i) {
          cheap_subsampling_ci(est[i], boot_est[, i], size_cb, n_val, alpha)
        }))
      },
      error = function(e) {
        stop(
          "Computation of confidence intervals failed with error: ",
          conditionMessage(e),
          ". Does your function return a vector of coefficients?"
        )
      }
    )
  }

  res$parameter <- rownames(res)
  rownames(res) <- NULL
  res <- list(
    res = res,
    b = b,
    size = size,
    type = type,
    alpha = alpha,
    n = n_val
  )
  if (keep_boot_estimates) {
    res <- c(res, list(boot_estimates = boot_est))
  }
  class(res) <- "cheap_bootstrap"
  res
}

##' Summary method for cheap_subsampling objects
##'
##' @title Summary method for cheap_subsampling objects
##' @param object An object of class "cheap_subsampling"
##' @param print Logical. Should the summary be printed? Defaults to TRUE.
##' @param ... Not applicable.
##' @return Summary table with the point estimates and confidence intervals.
##' @export
summary.cheap_bootstrap <- function(object, print = TRUE, ...) {
  if (print[[1]] == TRUE) {
    print(object, ...)
  }
  invisible(object$res)
}

##' Print method for cheap_bootstrap objects
##'
##' @title Print method for cheap_bootstrap objects
##' @param x An object of class "cheap_bootstrap"
##' @param ... Not applicable.
##' Prints the point estimates and confidence intervals.
##' @export
print.cheap_bootstrap <- function(x, ...) {
  if (x$type == "subsampling") {
    cat(
      paste0(
        "Cheap subsampling results for subsample size m = ",
        x$size,
        " and ",
        x$b,
        " bootstrap samples\n"
      )
    )
  } else {
    cat(
      paste0(
        "Cheap (non-parametric) bootstrap results for ",
        x$b,
        " bootstrap samples\n"
      )
    )
  }
  print(x$res, ...)
  invisible(x)
}

##' Get all bootstrap confidence intervals
##'
##' @title Get all bootstrap confidence intervals
##' @param x An object of class "cheap_bootstrap"
##' @param ... Not applicable.
##' Plots the point estimates and confidence intervals
##' as a function of the number of bootstrap samples.
##' @examples
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' cs <- cheap_bootstrap(x, b = 10, data = anorexia)
##' get_bootstrap_cis(cs)
##' @export get_bootstrap_cis.cheap_bootstrap
##' @export

##' @export
get_bootstrap_cis <- function(x, ...) {
  UseMethod("get_bootstrap_cis")
}

##' @rdname get_bootstrap_cis
##' @export get_bootstrap_cis.cheap_bootstrap
##' @export
get_bootstrap_cis.cheap_bootstrap <- function(x, ...) {
  est <- x$res$estimate
  if (is.null(x$boot_estimates)) {
    stop("No bootstrap estimates found in x")
  }
  res_b <- list()
  for (b_cur in seq_len(x$b)) {
    ## apply cheap_subsampling_ci for each row of boot_est, est

    tryCatch(
      {
        res <- do.call("rbind", lapply(seq_len(length(est)), function(i) {
          cheap_subsampling_ci(
            est[i],
            x$boot_est[seq_len(b_cur), i],
            x$size,
            x$n,
            x$alpha
          )
        }))
      },
      error = function(e) {
        stop(
          "Computation of confidence intervals failed with error: ",
          conditionMessage(e),
          " for b = ",
          b_cur
        )
      }
    )
    res$parameter <- x$res$parameter
    rownames(res) <- NULL
    res$b <- b_cur
    res_b[[b_cur]] <- res
  }
  do.call(rbind, res_b)
}

##' Plot method for cheap_bootstrap objects
##'
##' @title Plot method for cheap_bootstrap objects
##' @param x An object of class "cheap_bootstrap"
##' @param ... Not applicable.
##' Plots the point estimates and confidence intervals
##' as a function of the number of bootstrap samples.
##' @export
##' @examples
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' cs <- cheap_bootstrap(x, b = 10, data = anorexia)
##' plot(cs)
##'
plot.cheap_bootstrap <- function(x, ...) {
  b <- estimate <- cheap_lower <- cheap_upper <- NULL
  res_b <- get_bootstrap_cis(x, ...)
  ggplot2::ggplot(data = res_b, ggplot2::aes(x = b, y = estimate)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(
      alpha = 0.2,
      ggplot2::aes(
        ymin = cheap_lower,
        ymax = cheap_upper
      )
    ) +
    ggplot2::facet_wrap(~parameter, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Number of bootstrap samples") +
    ggplot2::ylab(paste0(
      "Cheap ",
      ifelse(x$type == "subsampling", "subsampling", "bootstrap"),
      " confidence intervals"
    ))
}
