# Constant predictors cannot help fit a nuisance model and may cause a learner
# to fail. Assess them in the actual fitting sample, including within strata.
continuous_varying_covariates <- function(data, covariates) {
    varies <- vapply(data[, covariates, with = FALSE],
                     function(x) data.table::uniqueN(x) > 1L, logical(1))
    covariates[varies]
}
