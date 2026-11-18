# Run from the rtmle repository root:
#   Rscript sandbox/test_continuous_time_sandbox.R
# Or in R: source("sandbox/test_continuous_time_sandbox.R")
#
# This is deliberately outside tests/testthat and excluded by .Rbuildignore.
# It loads BOTH working checkouts, fits each package to copies of the same
# data, and prints the two values and their differences. When sourced, the
# fitted objects are available in fits[[scenario]] for further inspection.
# Set CONTICEIPCW_PATH
# if the contICEIPCW checkout is not the sibling ../contICEIPCW.
# Source scenarios: contICEIPCW/tests/testthat/test_continuous_time_function.R.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
root <- if (length(file_arg)) {
    normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
} else normalizePath(getwd())
source_root <- normalizePath(Sys.getenv("CONTICEIPCW_PATH", file.path(root, "..", "contICEIPCW")), mustWork = TRUE)
if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload to run this sandbox.")
pkgload::load_all(source_root, export_all = FALSE, quiet = TRUE)
pkgload::load_all(root, export_all = FALSE, quiet = TRUE)

# The source simulator requires Rfast, which may not load on all systems.
# The migrated simulator is used only to create a common input; each analysis
# below receives an independent copy of the exact same data.
make_data <- function(uncensored, competing) {
    set.seed(34)
    rtmle::simulate_continuous_time_data(1000, uncensored = uncensored,
                                          no_competing_events = !competing)
}
copy_data <- function(d) lapply(d, data.table::copy)

source_fit <- function(d, horizons = 720, marginal = TRUE, source_learner = "oipcw_expit",
                       exclude = NULL, stratify = FALSE, intervene_all = FALSE,
                       censoring_formula = NULL, reduce_time = FALSE) {
    prepared <- contICEIPCW::prepare_data(copy_data(d), time_horizons = horizons,
        time_covariates = c("A", "L"), baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = marginal)
    scores <- contICEIPCW::propensity_scores(prepared, model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph", exclude_latest_covariate = exclude,
        marginal_censoring_formula = censoring_formula, reduce_colinearity_time = reduce_time)
    contICEIPCW::debias_ice_ipcw(scores, model_pseudo_outcome = source_learner,
        model_hazard = "learn_coxph", conservative = TRUE, return_ic = TRUE,
        stratify_sequential_regression = stratify,
        intervene_all_sequential_regression = intervene_all, reduce_colinearity_time = reduce_time)
}

rtmle_fit <- function(d, horizons = 720, marginal = TRUE, exclude = NULL,
                      stratify = FALSE, intervene_all = FALSE, censoring_formula = NULL,
                      learner = "learn_oipcw", scale_pseudo_outcome = NULL) {
    d <- copy_data(d)
    x <- rtmle::rtmle_init(c(0, horizons), name_id = "id", name_time = "time",
        name_outcome = "Y", name_competing = "D", name_censoring = "C", continuous_time = TRUE)
    x <- rtmle::add_baseline_data(x, d$baseline_data)
    x <- rtmle::add_continuous_data(x, event_data = d$timevarying_data,
        time_covariates = "L", treatment_events = "A", covariate_events = "L")
    x <- rtmle::prepare_rtmle_data(x)
    x <- rtmle::regime(x, "Always_A", intervention = 1, treatment_variables = "A")
    x <- rtmle::target(x, "risk", regimes = "Always_A")
    x <- rtmle::model_formula(x, verbose = FALSE)
    rtmle::run_rtmle(x, learner = learner, propensity_learner = "learn_glm",
        censoring_learner = list(name = "coxph", fun = "learn_coxph",
                                 marginal_censoring = marginal,
                                 marginal_censoring_formula = censoring_formula),
        time_horizon = seq_along(horizons),
        estimator = "one-step", exclude_latest_covariate = exclude,
        stratify_sequential_regression = stratify,
        intervene_all_sequential_regression = intervene_all,
        scale_pseudo_outcome = scale_pseudo_outcome)
}

# Unpenalized Gaussian regression is included as an additional direct check.
# It is not equivalent to the source's penalized LM test.
gaussian_learner <- function(character_formula, data, intervened_data, ...) {
    fit <- stats::glm(stats::as.formula(character_formula), data = data, family = stats::gaussian())
    list(predicted_values = stats::predict(fit, newdata = intervened_data, type = "response"))
}

# Each entry is a separate pair of fits. The first six cover the main
# censoring and competing-event settings in the source test file.
cases <- list(
    list(name = "uncensored", uncensored = TRUE, competing = FALSE),
    list(name = "uncensored, competing events", uncensored = TRUE, competing = TRUE),
    list(name = "conditional censoring", uncensored = FALSE, competing = FALSE, marginal = FALSE),
    list(name = "conditional censoring, competing events", uncensored = FALSE, competing = TRUE, marginal = FALSE),
    list(name = "marginal censoring", uncensored = FALSE, competing = FALSE),
    list(name = "marginal censoring, competing events", uncensored = FALSE, competing = TRUE),
    list(name = "intervene on all regression histories", uncensored = FALSE, competing = FALSE, marginal = FALSE, intervene_all = TRUE),
    list(name = "stratified sequential regression", uncensored = FALSE, competing = FALSE, marginal = FALSE, stratify = TRUE),
    list(name = "exclude latest time from propensity", uncensored = FALSE, competing = TRUE, exclude = "time"),
    list(name = "two horizons and exclude latest time", uncensored = FALSE, competing = TRUE, exclude = "time", horizons = c(40, 720)),
    list(name = "two horizons with conditional censoring", uncensored = FALSE, competing = TRUE,
         exclude = "time", horizons = c(40, 720), marginal = FALSE),
    list(name = "custom marginal censoring formula", uncensored = FALSE, competing = FALSE, censoring_formula = "~ age + A_0 + L_0"),
    list(name = "source time reparameterization (this fixture)", uncensored = FALSE, competing = TRUE, reduce_time = TRUE),
    list(name = "unpenalized Gaussian learner (additional)", uncensored = FALSE, competing = TRUE,
         source_learner = "lm", learner = gaussian_learner, scale_pseudo_outcome = FALSE),
    list(name = "covariate changes at treatment events", uncensored = FALSE, competing = TRUE,
         joint_covariate_events = TRUE)
)

# Absolute tolerance for these deterministic fits. Intervals deliberately use
# each package's own convention: source uses 1.96; rtmle uses qnorm(.975).
tolerance <- 1e-7
failures <- character()
fits <- list()
for (case in cases) {
    name <- case$name
    cat("\n---", name, "---\n")
    tryCatch({
        d <- make_data(case$uncensored, case$competing)
        if (isTRUE(case$joint_covariate_events)) {
            # Each treatment event changes L by a subject-dependent amount;
            # retain that change at all subsequent rows in the same history.
            # Terminal rows also receive distinct recorded covariate values.
            d$timevarying_data[, L := L + cumsum((event == "A") * (1 + id %% 3)) +
                                  (event %in% c("C", "Y", "D", "tauend")) * (2 + id %% 5), by = id]
        }
        shared <- case[intersect(names(case), c("horizons", "marginal", "exclude", "stratify", "intervene_all", "censoring_formula"))]
        source_args <- c(list(d = d), shared, case[intersect(names(case), c("source_learner", "reduce_time"))])
        target_args <- c(list(d = d), shared, case[intersect(names(case), c("learner", "scale_pseudo_outcome"))])
        source <- do.call(source_fit, source_args)
        target <- do.call(rtmle_fit, target_args)
        fits[[name]] <- list(contICEIPCW = source, rtmle = target)
        if (name %in% c("two horizons and exclude latest time",
                        "two horizons with conditional censoring")) {
            # Compare ordinary and pooled interval predictions as well as the
            # final estimates, aligning every subject before comparing.
            marginal <- if (is.null(case$marginal)) TRUE else case$marginal
            source_prepared <- contICEIPCW::prepare_data(copy_data(d),
                time_horizons = case$horizons, time_covariates = c("A", "L"),
                baseline_covariates = c("age", "A_0", "L_0"), marginal_censoring = marginal)
            source_scores <- contICEIPCW::propensity_scores(source_prepared,
                model_treatment = "learn_glm_logistic", model_hazard = "learn_coxph",
                exclude_latest_covariate = case$exclude)
            source_info <- source_prepared$info[order(time_horizon)]
            target_info <- target$continuous_prepared$info[order(time_horizon)]
            if (!identical(as.numeric(source_info$time_horizon), as.numeric(target_info$time_horizon)) ||
                !identical(as.integer(source_info$last_event), as.integer(target_info$last_event)) ||
                length(unique(source_info$last_event)) < 2L)
                stop("Multiple-horizon event cutoffs differ or are not distinct")
            source_nuisance <- source_scores$data
            target_nuisance <- target$continuous_nuisance$Always_A$data
            index <- match(source_nuisance$id, target_nuisance$id)
            if (anyNA(index) || anyDuplicated(index) || length(index) != nrow(target_nuisance))
                stop("Cannot align censoring predictions by subject")
            columns <- c(paste0("survival_censoring_", seq_len(max(source_info$last_event))),
                         paste0("survival_censoring_pooled_", unique(source_info$last_event)))
            differences <- lapply(columns, function(column) {
                if (!column %in% names(source_nuisance) || !column %in% names(target_nuisance))
                    stop("Missing censoring prediction column: ", column)
                a <- as.numeric(source_nuisance[[column]])
                b <- as.numeric(target_nuisance[[column]][index])
                if (!identical(is.na(a), is.na(b))) stop("Censoring prediction rows differ: ", column)
                observed <- !is.na(a)
                difference <- if (any(observed)) max(abs(a[observed] - b[observed])) else 0
                if (!is.finite(difference) || difference > tolerance)
                    stop("Censoring predictions differ: ", column)
                data.frame(column = column, subjects = sum(observed), max_difference = difference)
            })
            cat("Event cutoffs by horizon:\n")
            print(source_info[, c("time_horizon", "last_event"), with = FALSE])
            cat(if (marginal) "Marginal" else "Conditional", "censoring prediction differences:\n")
            print(do.call(rbind, differences), row.names = FALSE, digits = 6)
        }
        if (isTRUE(case$joint_covariate_events)) {
            # Compare all retained event marks and covariate values, not just
            # estimates, which could agree despite a missed update.
            source_history <- contICEIPCW::prepare_data(copy_data(d), time_horizons = 720,
                time_covariates = c("A", "L"), baseline_covariates = c("age", "A_0", "L_0"),
                marginal_censoring = TRUE)$wide_data
            target_history <- target$continuous_prepared$wide_data
            index <- match(source_history$id, target_history$id)
            if (anyNA(index) || anyDuplicated(index)) stop("Cannot align prepared histories by ID")
            checked <- treatment_changed <- terminal_changed <- 0L
            events <- grep("^event_[1-9][0-9]*$", names(source_history), value = TRUE)
            for (event_column in events) {
                k <- sub("^event_", "", event_column)
                L_column <- paste0("L_", k)
                previous_L <- paste0("L_", as.integer(k) - 1L)
                marks <- as.character(source_history[[event_column]])
                rows <- which(!is.na(marks))
                if (!length(rows)) next
                if (!all(c(event_column, L_column) %in% names(target_history)) ||
                    !identical(as.character(target_history[[event_column]][index[rows]]), marks[rows]) ||
                    !isTRUE(all.equal(as.numeric(target_history[[L_column]][index[rows]]),
                                      as.numeric(source_history[[L_column]][rows]), tolerance = tolerance)))
                    stop("Event marks or covariate histories differ")
                checked <- checked + length(rows)
                changed <- source_history[[L_column]][rows] != source_history[[previous_L]][rows]
                treatment_changed <- treatment_changed + sum(changed & marks[rows] == "A")
                terminal_changed <- terminal_changed + sum(changed & marks[rows] %in% c("C", "Y", "D", "tauend"))
            }
            if (!treatment_changed || !terminal_changed)
                stop("Fixture must contain covariate changes at treatment and terminal events")
            cat("Event records checked:", checked, "; treatment events with a changed covariate:",
                treatment_changed, "; terminal events with a changed covariate:", terminal_changed, "\n")
        }
        a <- as.data.frame(source$results)
        b <- as.data.frame(target$estimate$Main_analysis)
        a <- a[order(a$time_horizon), ]
        b <- b[order(b$Time_horizon), ]
        target_horizons <- target$time_grid_scale[b$Time_horizon + 1L]
        if (!identical(as.numeric(a$time_horizon), as.numeric(target_horizons))) stop("Horizon mismatch")
        source_values <- cbind(estimate = a$estimate, se = a$se, ice_ipcw = a$ice_ipcw_estimate, ipw = a$ipw)
        rtmle_values <- cbind(estimate = b$Estimate, se = b$Standard_error, ice_ipcw = b$ICE_IPCW, ipw = b$IPW)
        tab <- data.frame(horizon = rep(a$time_horizon, each = ncol(source_values)),
            quantity = rep(colnames(source_values), times = nrow(source_values)),
            contICEIPCW = as.vector(t(source_values)), rtmle = as.vector(t(rtmle_values)))
        tab$difference <- tab$rtmle - tab$contICEIPCW
        print(tab, row.names = FALSE, digits = 10)
        if (any(!is.finite(tab$difference)) || any(abs(tab$difference) > tolerance)) stop("Estimate/SE/ICE/IPW comparison failed")
        # Align individual influence curves by subject ID, not by assumed row order.
        ids <- as.character(target$continuous_prepared$wide_data$id)
        for (h in a$time_horizon) {
            source_ic <- source$ic$ic[[match(h, source$ic$time_horizons)]]
            node <- match(h, target$time_grid_scale) - 1L
            target_ic <- target$IC$risk$Always_A[[paste0("time_horizon_", node)]]
            source_ids <- if (is.matrix(source_ic)) rownames(source_ic) else names(source_ic)
            index <- match(ids, source_ids)
            if (length(target_ic) != length(ids) || anyNA(index) || anyDuplicated(ids))
                stop("Cannot align individual influence curves by ID")
            ic_difference <- max(abs(as.numeric(source_ic)[index] - target_ic))
            cat("Horizon", h, "maximum individual influence-curve difference:", format(ic_difference, digits = 6), "\n")
            if (!is.finite(ic_difference) || ic_difference > tolerance) stop("Individual influence curves differ")
        }
        # Validate intervals separately because their critical values differ.
        if (any(abs(a$lower - (a$estimate - 1.96 * a$se)) > tolerance) ||
            any(abs(a$upper - (a$estimate + 1.96 * a$se)) > tolerance) ||
            any(abs(b$Lower - (b$Estimate - stats::qnorm(.975) * b$Standard_error)) > tolerance) ||
            any(abs(b$Upper - (b$Estimate + stats::qnorm(.975) * b$Standard_error)) > tolerance)) stop("Confidence interval convention check failed")
        cat("PASS (absolute tolerance", tolerance, "; confidence intervals checked using each package's critical value)\n")
    }, error = function(e) {
        failures <<- c(failures, paste0(name, ": ", conditionMessage(e)))
        cat("FAIL:", conditionMessage(e), "\n")
    })
}
cat("\nNot asserted equal: source TMLE drops the final baseline update; source treatment-0 weights and propensity bounding have known bugs; conservative = FALSE is not implemented in rtmle.\n")
cat("Not claimed equivalent: source penalized, probit, or nonlinear learners and source-specific coefficient output. The time reparameterization comparison applies to this unpenalized fixture only.\n")
cat("The source data_large tests use different lag settings for propensity and outcome; rtmle currently exposes one shared lag, so those exact specifications are not compared.\n")
if (length(failures)) stop(length(failures), " comparison(s) failed:\n", paste(failures, collapse = "\n"))
cat("\nAll", length(cases), "side-by-side scenarios passed.\n")
