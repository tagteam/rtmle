# Censoring model structure belongs to its learner specification, not data
# preparation. Remove these model options before invoking the learner itself.
continuous_censoring_specification <- function(learner) {
    if (is.null(learner)) return(NULL)
    parsed <- continuous_nuisance_learner(learner)
    if (identical(parsed$name, "superlearn"))
        stop("The censoring learner must be a single survival learner; censoring ensembles are not supported.")
    time_varying <- parsed$args[["time_varying", exact = TRUE]]
    if (is.null(time_varying)) time_varying <- FALSE
    if (!is.logical(time_varying) || length(time_varying) != 1L || is.na(time_varying))
        stop("censoring_learner time_varying must be TRUE or FALSE.")
    marginal <- parsed$args[["marginal_censoring", exact = TRUE]]
    if (is.null(marginal)) marginal <- !time_varying
    if (!is.logical(marginal) || length(marginal) != 1L || is.na(marginal))
        stop("censoring_learner marginal_censoring must be TRUE or FALSE.")
    if (time_varying && marginal)
        stop("time_varying = TRUE cannot be combined with marginal_censoring = TRUE.")
    formula <- parsed$args[["marginal_censoring_formula", exact = TRUE]]
    if (!is.null(formula)) {
        if (!marginal) stop("marginal_censoring_formula requires marginal_censoring = TRUE.")
        if (inherits(formula, "formula")) formula <- paste(deparse(formula), collapse = " ")
        if (!is.character(formula) || length(formula) != 1L || is.na(formula) ||
            !startsWith(trimws(formula), "~"))
            stop("marginal_censoring_formula must be a one-sided formula.")
    }
    parsed$args[["marginal_censoring"]] <- NULL
    parsed$args[["marginal_censoring_formula"]] <- NULL
    current_formula <- parsed$args[["censoring_formula", exact = TRUE]]
    if (!is.null(current_formula)) {
        if (!time_varying) stop("censoring_formula requires time_varying = TRUE.")
        if (inherits(current_formula, "formula")) current_formula <- paste(deparse(current_formula), collapse = " ")
        if (!is.character(current_formula) || length(current_formula) != 1L || is.na(current_formula) ||
            !startsWith(trimws(current_formula), "~"))
            stop("censoring_formula must be a one-sided formula.")
    }
    parsed$args[["time_varying"]] <- NULL
    parsed$args[["censoring_formula"]] <- NULL
    list(learner = parsed, marginal = marginal, formula = formula,
         time_varying = time_varying, current_formula = current_formula)
}
