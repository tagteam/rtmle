continuous_cox_survival <- function(fit, data, time_variable = "time",
                                    time_ref = NULL, left_limit = TRUE) {
    endpoint <- data[[time_variable]]
    start <- if (is.null(time_ref)) rep(0, nrow(data)) else data[[time_ref]]
    if (!is.numeric(endpoint) || length(endpoint) != nrow(data) ||
        !is.numeric(start) || length(start) != nrow(data) ||
        any(!is.finite(endpoint)) || any(!is.finite(start)) ||
        any(start < 0 | endpoint < start)) {
        stop("Censoring interval start and endpoint must be finite, nonnegative, and ordered.")
    }
    if (is.null(fit$fit)) return(rep(1, nrow(data)))
    hazard <- fit$baseline_hazard
    # findInterval gives right-continuous values; exact endpoint matches must
    # step back once to evaluate the left limit. Each Cox stratum has its own
    # cumulative hazard and its own ordered time scale.
    if ("strata" %in% names(hazard)) {
        terms <- stats::delete.response(stats::terms(fit$fit))
        special <- attr(terms, "specials")$strata
        variables <- rownames(attr(terms, "factors"))[special]
        if (length(variables) != 1L)
            stop("Use a single strata() term, combining variables inside strata() if needed.")
        frame <- stats::model.frame(terms, data = data, xlev = fit$fit$xlevels,
                                    na.action = stats::na.fail)
        groups <- as.character(frame[[variables]])
        hazard_groups <- as.character(hazard$strata)
    } else {
        groups <- rep("all", nrow(data))
        hazard_groups <- rep("all", nrow(hazard))
    }
    increment <- numeric(nrow(data))
    for (group in unique(groups)) {
        rows <- which(groups == group)
        group_hazard <- hazard[hazard_groups == group, , drop = FALSE]
        if (!nrow(group_hazard)) stop("Censoring prediction includes an unknown Cox stratum.")
        end_index <- findInterval(endpoint[rows], group_hazard$time)
        if (left_limit) end_index <- end_index - as.integer(endpoint[rows] %in% group_hazard$time)
        start_index <- findInterval(start[rows], group_hazard$time)
        h <- c(0, group_hazard$hazard)
        increment[rows] <- pmax(0, h[end_index + 1L] - h[start_index + 1L])
    }
    # Under proportional hazards, S(end- | start) is exp(-exp(lp) times
    # the cumulative baseline hazard increment over the interval).
    lp <- stats::predict(fit$fit, newdata = data, type = "lp", reference = "zero")
    predictions <- exp(-increment * exp(lp))
    continuous_nuisance_predictions(list(predicted_values = predictions), nrow(data), "Censoring", TRUE)
}
