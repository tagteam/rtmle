#' Perform the Cheap Subsampling Procedure for LTMLE
#'
#' This function applies a bootstrap resampling technique using the `CheapSubsampling` package.
#' It modifies the provided dataset `x` by resampling observations and then applies LTMLE.
#'
#' @param x An object contaiining a prepared `rtmle` object
#' @param cheap_bootstrap_arguments A list specifying bootstrap parameters given to the `cheap_bootstrap` function.
#'   - `b` (integer): Number of bootstrap samples (default: 25).
#'   - `type` (character): Type of resampling method, e.g., "subsampling".
#'   - `size` (integer, optional): Size of each subsample
#' @param rtmle_arguments A list of arguments passed to the `run_rtmle` function.
#' @param id_name A character string specifying the name of the identifier column used for personal ids.
#'
#' @examples
#'
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
#' cb <- cheap_bootstrap_rtmle(
#'  x = x,
#'  rtmle_arguments = list(learner = "learn_glm",time_horizon = tau),
#'  id_name = "id",
#'  cheap_bootstrap_arguments = list(b = 3, type = "subsampling", size = NULL)
#' )
#' x <- run_rtmle(x,learner = "learn_glm",time_horizon = tau)
#' summary(x)
#' cb
#'
#' @export
cheap_bootstrap_rtmle <- function(x,
                                  cheap_bootstrap_arguments = list(b=25, type = "subsampling", size = NULL),
                                  rtmle_arguments,
                                  id_name) {
  ## Recursively apply rbind
  recursive_rbind <- function(x) if (inherits(x, "list")) do.call("rbind", lapply(x, recursive_rbind)) else x

  fun <- function(data) {
    ## Assume that the ids across x$followup and x$prepared_data are ordered the same way
    matched_ids <- match(data[[id_name]], x$prepared_data[[id_name]])

    ## Subset data according to the resampled data]
    x$prepared_data <- x$prepared_data[matched_ids]
    x$followup <- x$followup[matched_ids]

    ## Run rtmle and suppress messages
    suppressMessages(res <- recursive_rbind(summary(do.call(
      "run_rtmle",
      append(list(x = x), rtmle_arguments)
    ))))

    ## For correct usage of CheapSubsampling package, no numerics are allowed in paramater_column_names
    ## Thus Time_horizon should be changed to, say, a factor
    res$Time_horizon <- as.factor(res$Time_horizon)
    res[, c("Target", "Protocol", "Target_parameter", "Time_horizon", "Estimator","Estimate")]
  }

  do.call(
    "cheap_bootstrap",
    append(
      list(
        fun = fun,
        estimate_column_name = "Estimate",
        parameter_column_names = c("Target", "Protocol", "Target_parameter", "Time_horizon", "Estimator"),
        data = x$prepared_data
      ),
      cheap_bootstrap_arguments
    )
  )
}
## Retrieve cheap subsampling confidence interval
confidence_interval <- function(est,
                                est_boot,
                                size,
                                sample_size,
                                alpha,
                                type = "subsampling") {
  b_max <- length(est_boot)
  res <- list()
  for (b in seq_len(b_max)) {
    b_est <- est_boot[seq_len(b)]
    if (type == "subsampling") {
      sub_factor <- sqrt((size) / (sample_size - size))
    } else {
      sub_factor <- 1
    }
    s_val <- sqrt(mean((est - b_est)^2))
    tq <- stats::qt(1 - alpha / 2, df = b)
    res[[b]] <- data.frame(
      estimate = est,
      cheap_lower = est - tq * sub_factor * s_val,
      cheap_upper = est + tq * sub_factor * s_val,
      b = b
    )
  }
  do.call(rbind, res)
}

##' Method implementing the cheap subsampling method for confidence intervals
##'
##' Given a model object or a function that returns a vector of coefficients
##' and a data set,
##' this function computes confidence intervals using
##' the cheap subsampling method.
##'
##' @title Cheap subsampling
##' @param fun A function that returns a data.frame
## 'with the point estimates given
##' in estimate_column_name and the corresponding parameter names in parameter_column_names.
##' @param data Data set to be used for the computation.
##' @param parameter_column_names Character. Column names of the parameters.
##' @param estimate_column_name Character. Column name of the point estimates.
##' @param b Number of bootstrap samples.
##' @param size Subsample size. Defaults to 0.632 * nrow(data).
##' @param alpha Significance level. Defaults to 0.05.
##' @param type Character. Type of bootstrap method.
##' Can be either "subsampling" or "non_parametric".
##' Defaults to "subsampling".
##' @param parallelize Logical. Should the computation be parallelized?
##' @param cores Number of cores to be used for parallel computation.
##' @param keep_estimates Logical. Should the bootstrapped estimates be kept?
##' @return An object of class "cheap_bootstrap" containing
##' the point estimates and confidence intervals.
##' @export
##' @examples
##' \dontrun{
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' library(broom)
##' library(dplyr)
##' fun <- function(data) {
##'   lm(formula = Postwt ~ Prewt + Treat + offset(Prewt), data = data) %>%
##'   tidy()
##' }
##' cs <- cheap_bootstrap(fun = fun,
##'                       b = 20,
##'                       data = anorexia,
##'                       estimate_column_name = "estimate",
##'                       parameter_column_names = "term")
##' cs
##' summary(cs)
##' plot(cs)
##' }
##' ## example with a function call and parallel computation
##' \dontrun{
##' set.seed(123)
##' ## note the function needs to load the packages needed
##' ## for the computation if parallelized
##' cs <- cheap_bootstrap(fun = fun,
##'                       b = 20,
##'                       data = anorexia,
##'                       estimate_column_name = "estimate",
##'                       parameter_column_names = "term",
##'                       parallelize = TRUE,
##'                       cores = 2)
##' }
cheap_bootstrap <- function(fun,
                            data,
                            size = NULL,
                            type = "subsampling",
                            b = 20,
                            estimate_column_name = "estimate",
                            parameter_column_names = "parameter",
                            alpha = 0.05,
                            parallelize = FALSE,
                            cores = 1,
                            keep_estimates = TRUE) {
  estimate <- b_est <- NULL
  ## Check if the data is a data.frame
  if (!inherits(data, "data.frame")) {
    stop("data must inherit data.frame")
  }
  sample_size <- nrow(data)

  ## If parallelize, check that namespace is loaded
  if (parallelize) {
    requireNamespace("parallel")
  }

  ## Check that all relevant arguments are of length 1
  arguments_check <- c(
    "b",
    "alpha",
    "type",
    "parallelize",
    "cores",
    "keep_estimates"
  )
  for (variable_name in arguments_check) {
    if (length(get(eval(variable_name))) != 1) {
      stop(paste(variable_name, "must be of length 1"))
    }
  }

  if (is.null(size) && type == "subsampling") {
    size <- round(0.632 * sample_size)
  } else if (type == "non_parametric") {
    size <- sample_size
  }

  est <- tryCatch(
    {
      fun(data)
    },
    error = function(e) {
      stop(
        "function fun failed with error on the original dataset: ",
        conditionMessage(e)
      )
    }
  )
  ## Check that est inherits data.frame
  if (!inherits(est, "data.frame")) {
    stop("function fun/derived from fun must return a data.frame")
  }

  est <- data.table::as.data.table(est)
  original_value <- est

  ## Check that estimate_column_name and parameter_column_names
  ## are in the names of columns of est
  if (!(estimate_column_name %in% colnames(est)) ||
      !(all(parameter_column_names %in% colnames(est)))) {
    stop("estimate_column_name and parameter_column_names
         must be in the names of columns of est")
  }

  ## Check that that the column specified by estimate_column_name is numeric
  if (!all(sapply(est[[estimate_column_name]], is.numeric))) {
    stop("The column specified by estimate_column_name must be numeric")
  }

  ## Check that the column specified by parameter_column_names is character or factor
  if (sum(est[, sapply(.SD, function(x) !is.character(x) && !is.factor(x)),
              .SDcols = parameter_column_names
  ]) > 0) {
    stop("The column specified by parameter_column_names must be character or factor")
  }
  est <- est[, c(parameter_column_names, estimate_column_name), with = FALSE]

  boot_fun <- function(i) {
    set.seed(seeds[i])
    tryCatch(
      {
        x <-
          fun(data[sample(1:sample_size, size, replace = (type == "non_parametric")), ,
                   drop = FALSE
          ])
        x$b <- i
        x
      },
      error = function(e) {
        stop(
          "Bootstrap computation failed with error: ",
          conditionMessage(e),
          " for iteration ", i, " with the seed ", seeds[i]
        )
      }
    )
  }

  ## sample b seeds
  seeds <- sample.int(1e+09, b)

  if (parallelize) {
    results <- pbmcapply::pbmclapply(
      X = seq_len(b),
      FUN = boot_fun,
      mc.cores = cores
    )
  } else {
    results <- lapply(seq_len(b), boot_fun)
  }
  tryCatch(
    {
      boot_est <- data.table::as.data.table(do.call("rbind", results))[, c(parameter_column_names, estimate_column_name, "b"), with = FALSE]
    },
    error = function(e) {
      stop("Could not bind the results.
            Did you return a data.frame in the function fun?")
    }
  )

  ## Rename the column by estimate_column_name to b_est
  ## and rename the column from est to estimate
  data.table::setnames(boot_est, estimate_column_name, "b_est")
  data.table::setnames(est, estimate_column_name, "estimate")

  ## Merge boot_est with by parameter_column_names
  boot_est <- merge(est,
                    boot_est,
                    by = parameter_column_names,
                    all.x = TRUE
  )

  boot_est[, c("size", "sample_size", "alpha", "type") := list(size, sample_size, alpha, type)]

  tryCatch(
    {
      ## Apply confidence_interval for each column specified by parameter_column_names
      res <- boot_est[, confidence_interval(
        est = estimate[1],
        est_boot = b_est,
        size = size[1],
        sample_size = sample_size[1],
        alpha = alpha[1],
        type = type[1]
      ),
      by = parameter_column_names
      ]
    },
    error = function(e) {
      stop(
        "Computation of confidence intervals failed with error: ",
        conditionMessage(e),
        ". Does your function return a data.frame?"
      )
    }
  )

  res <- list(
    res = res,
    b = b,
    size = size,
    type = type,
    alpha = alpha,
    n = sample_size,
    parameter_column_names = parameter_column_names,
    original_value = original_value,
    seeds = seeds
  )
  if (keep_estimates) {
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

##' Plot method for cheap_bootstrap objects
##'
##' @title Plot method for cheap_bootstrap objects
##' @param x An object of class "cheap_bootstrap"
##' @param alternate_confidence_interval Numeric.
##' Additional confidence interval to be added to the plot.
##' Must be a data.frame consisting of the columns specified by parameter_column_names
##' and two additional columns of type numeric.
##' for each parameter.
##' @param ... Not applicable.
##' Plots the point estimates and confidence intervals
##' as a function of the number of bootstrap samples.
##' @export
plot.cheap_bootstrap <- function(x, alternate_confidence_interval = NULL, ...) {
  form_par_wrap <- stats::as.formula(paste0("~", paste0(x$parameter_column_names, collapse = "+")))
  b <- estimate <- cheap_lower <- cheap_upper <- lower <- upper <- NULL
  if (!is.null(alternate_confidence_interval)) {
    if (!inherits(alternate_confidence_interval, "data.frame")) {
      stop("alternate_confidence_interval must be a data.frame")
    }
    alternate_confidence_interval <- data.table::as.data.table(alternate_confidence_interval)

    names_extra <- colnames(alternate_confidence_interval)
    names_lower_upper <- setdiff(names_extra, x$parameter_column_names)
    ## Check that x$parameter_column_names is in the names of alternate_confidence_interval and that
    ## two additional columns are available and of numeric type
    if (!(all(x$parameter_column_names %in% names_extra)) ||
        !(all(sapply(
          alternate_confidence_interval[, names_lower_upper, with = FALSE],
          is.numeric
        )))) {
      stop("alternate_confidence_interval must contain the columns specified by parameter_column_names
           and two additional columns of type numeric")
    }

    ## Merge alternate_confidence_interval with x$res by parameter_column_names
    x$res <- merge(x$res, alternate_confidence_interval, by = x$parameter_column_names, all.x = TRUE)

    data.table::setnames(x$res, names_lower_upper, c("lower", "upper"))
  }

  p <- ggplot2::ggplot(
    data = x$res,
    ggplot2::aes(x = b, y = estimate)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(
      alpha = 0.2,
      ggplot2::aes(
        ymin = cheap_lower,
        ymax = cheap_upper
      )
    ) +
    ggplot2::facet_wrap(form_par_wrap,
                        scales = "free_y",
                        labeller = ggplot2::label_both
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Number of bootstrap samples") +
    ggplot2::ggtitle(
      paste0(
        "Cheap ",
        ifelse(x$type == "subsampling", "Subsampling", "Bootstrap"),
        " confidence intervals compared to IC based intervals
        (subsample of size m = ",
        x$size,
        " out of n = ",
        x$n,
        " observations (",
        round(x$size / x$n * 100),
        "%)))"
      )
    ) +
    ggplot2::ylab("")
  if (!is.null(alternate_confidence_interval)) {
    p <- p +
      ggplot2::geom_hline(aes(yintercept = lower)) +
      ggplot2::geom_hline(aes(yintercept = upper))
  }
  p
}
