# Build ordinary additive nuisance models without making callers construct strings.
continuous_nuisance_formula <- function(outcome, covariates) {
    rhs <- if (length(covariates)) paste(sprintf("`%s`", covariates), collapse = " + ") else "1"
    paste(outcome, "~", rhs)
}
