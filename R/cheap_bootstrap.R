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
#' rtmle_arguments = list(learner = "learn_glm",time_horizon = 1:tau),
#'   id_name = "id"
#' )
#' x <- run_rtmle(x,learner = "learn_glm",time_horizon = 1:tau)
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
    suppressMessages(res <- recursive_rbind(unique(summary(do.call(
      "run_rtmle",
      append(list(x = x), rtmle_arguments)
    )))))

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
