# At a continuous treatment opportunity the current covariate and time marks
# are observed before assignment. The treatment at this node is not.
continuous_adherence_strata <- function(strata, node, available_names,
                                        baseline_covariates, time_covariates) {
    validate_adherence_model_strata(strata)
    if (!length(strata)) return(NULL)
    resolve <- function(reference) {
        if (node == 0L && reference %in% c("A", "A_0"))
            stop("`adherence_model_strata` must not use the baseline treatment being assigned.")
        if (reference %in% baseline_covariates) return(reference)
        if (reference == "A") return(paste0("A_", node - 1L))
        if (reference %in% c(setdiff(time_covariates, "A"), "time", "event"))
            return(paste0(reference, "_", node))
        if (grepl("_[0-9]+$", reference)) {
            index <- as.integer(sub("^.*_", "", reference))
            base <- sub("_[0-9]+$", "", reference)
            if ((base == "A" && index >= node) ||
                (base %in% c(time_covariates, "time", "event") && index > node))
                stop("`adherence_model_strata` must use pre-decision variables; ", reference,
                     " is not available before treatment at node ", node, ".")
        }
        reference
    }
    resolved <- unique(vapply(strata, resolve, character(1)))
    missing <- setdiff(resolved, available_names)
    if (length(missing)) stop("Unknown variable(s) supplied in `adherence_model_strata`: ",
                              paste(missing, collapse = ", "))
    resolved
}
