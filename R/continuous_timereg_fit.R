# Fit one additive or Cox-Aalen model. Predictors stay fixed within each row's
# observation interval; repeated rows for one person describe their updates.
continuous_timereg_fit <- function(character_formula, data, model, id_variable, args) {
    if (!requireNamespace("timereg", quietly = TRUE))
        stop("Install the 'timereg' package to use Aalen or Cox-Aalen censoring models.")
    if (!is.character(id_variable) || length(id_variable) != 1L || is.na(id_variable) || !nzchar(id_variable))
        stop("id_variable must name one subject identifier column.")
    data <- droplevels(as.data.frame(data))
    formula <- stats::as.formula(character_formula)
    calls <- continuous_formula_calls(formula)
    if (any(calls %in% c("|", "||")) || any(grepl("^frailty($|[.])", calls)))
        stop("Random effects are not supported by the censoring learners.")
    if ("cluster" %in% calls)
        stop("Clustered variance formulas are not supported by these censoring learners; use id_variable to identify repeated intervals.")
    formula_env <- new.env(parent = environment(formula))
    formula_env$prop <- timereg::prop
    formula_env$const <- timereg::const
    environment(formula) <- formula_env
    terms_text <- paste(deparse(formula), collapse = " ")
    if (model == "aalen" && grepl("prop\\s*\\(", terms_text))
        stop("Use learn_cox_aalen for proportional prop() terms.")
    if (model == "cox.aalen" && grepl("const\\s*\\(", terms_text))
        stop("Cox-Aalen uses prop() for proportional effects; const() is an Aalen term.")
    frame <- stats::model.frame(formula, data, na.action = stats::na.fail, drop.unused.levels = TRUE)
    response <- stats::model.response(frame)
    if (!inherits(response, "Surv") || !attr(response, "type") %in% c("right", "counting"))
        stop("Aalen and Cox-Aalen learners require a right-censored or start-stop Surv response.")
    status <- response[, ncol(response)]
    stop_time <- response[, ncol(response) - 1L]
    if (anyDuplicated(stop_time[status == 1]))
        stop("Aalen and Cox-Aalen censoring models currently require distinct censoring times across subjects: timereg otherwise breaks ties randomly. Use the Cox learner for tied censoring times.")
    counting <- identical(attr(response, "type"), "counting")
    if (counting && (!id_variable %in% names(data) || anyNA(data[[id_variable]])))
        stop("Start-stop data need a nonmissing subject identifier column named '", id_variable, "'.")
    if (any(stats::model.offset(frame) != 0) || (length(args$offsets) && any(args$offsets != 0)))
        stop("Nonzero offsets are not supported by the Aalen or Cox-Aalen censoring learners.")
    rhs <- stats::delete.response(stats::terms(frame))
    result <- list(fit = NULL, formula = formula, terms = rhs, model = model,
                   xlevels = stats::.getXlevels(stats::terms(frame), frame))
    if (!any(status == 1)) return(result)
    design_frame <- stats::model.frame(rhs, data, na.action = stats::na.fail, drop.unused.levels = TRUE)
    design <- timereg::read.design(design_frame, rhs, model = model)
    if (!design$px) stop("The additive part of the censoring model needs an intercept or a predictor.")
    # timereg otherwise silently skips events when the additive design cannot
    # be estimated from the people still at risk at that time.
    start_time <- if (counting) response[, 1L] else rep(0, nrow(data))
    positive_weight <- rep(TRUE, nrow(data))
    if (length(args$weights)) {
        if (!is.numeric(args$weights) || length(args$weights) != nrow(data) ||
            any(!is.finite(args$weights)) || any(args$weights < 0))
            stop("Censoring learner weights must be finite, nonnegative, and match the training rows.")
        positive_weight <- args$weights > 0
    }
    for (time in stop_time[status == 1]) {
        at_risk <- start_time < time & stop_time >= time & positive_weight
        if (qr(design$X[at_risk, , drop = FALSE])$rank < design$px)
            stop("The additive censoring predictors are not identifiable in the risk set at time ", time,
                 "; simplify censoring_formula. The model would otherwise silently skip this censoring event.")
    }
    # Without prop() terms the Cox-Aalen model is the pure Aalen model.
    if (model == "cox.aalen" && !design$pz) model <- "aalen"
    result$model <- model
    result$additive_names <- design$covnamesX
    result$fixed_names <- design$covnamesZ
    defaults <- list(formula = formula, data = as.data.frame(data), n.sim = 0,
                     robust = 0, resample.iid = 0)
    if (counting) defaults$id <- data[[id_variable]]
    fit_args <- utils::modifyList(defaults, args)
    if (!identical(fit_args$formula, formula) || !identical(fit_args$data, defaults$data))
        stop("Specify the model formula and data through the learner's named arguments.")
    fun <- if (model == "aalen") timereg::aalen else timereg::cox.aalen
    # cox.aalen consumes random numbers even with simulations disabled.
    # Preserve the caller's seed so it does not alter later learner fits.
    if (identical(fit_args$n.sim, 0) || identical(fit_args$n.sim, 0L)) {
        old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit(if (is.null(old_seed)) {
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                rm(".Random.seed", envir = .GlobalEnv)
        } else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    fit <- do.call(fun, fit_args)
    if (any(!is.finite(fit$cum)) || any(!is.finite(fit$gamma)))
        stop("The censoring model returned nonfinite coefficients; simplify the formula or check the risk sets.")
    result$fit <- fit
    result
}
