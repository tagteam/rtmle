# Integrate a fitted hazard with the row's covariates held fixed. The start is
# excluded; left_limit selects whether a jump at the endpoint is excluded too.
continuous_timereg_survival <- function(result, data, time_variable, time_ref,
                                        left_limit = TRUE) {
    if (!is.logical(left_limit) || length(left_limit) != 1L || is.na(left_limit))
        stop("left_limit must be TRUE or FALSE.")
    endpoint <- data[[time_variable]]
    start <- if (is.null(time_ref)) rep(0, nrow(data)) else data[[time_ref]]
    if (!is.numeric(endpoint) || length(endpoint) != nrow(data) ||
        !is.numeric(start) || length(start) != nrow(data) ||
        any(!is.finite(endpoint)) || any(!is.finite(start)) ||
        any(start < 0 | endpoint < start))
        stop("Censoring interval start and endpoint must be finite, nonnegative, and ordered.")
    if (is.null(result$fit)) return(rep(1, nrow(data)))
    data <- as.data.frame(data)
    for (variable in intersect(names(result$xlevels), names(data))) {
        data[[variable]] <- factor(as.character(data[[variable]]), levels = result$xlevels[[variable]])
        if (anyNA(data[[variable]])) stop("Unknown or missing factor level in censoring predictor '", variable, "'.")
    }
    rhs <- result$terms
    frame <- stats::model.frame(rhs, data, xlev = result$xlevels, na.action = stats::na.fail)
    design <- timereg::read.design(frame, rhs, model = result$model)
    if (!identical(design$covnamesX, result$additive_names) ||
        !identical(design$covnamesZ, result$fixed_names))
        stop("Censoring prediction design does not match the fitted model.")
    cumulative <- result$fit$cum
    times <- cumulative[, 1L]
    coefficients <- cumulative[, -1L, drop = FALSE]
    jumps <- rbind(coefficients[1L, , drop = FALSE],
                   coefficients[-1L, , drop = FALSE] - coefficients[-nrow(coefficients), , drop = FALSE])
    fixed <- if (design$pz) as.numeric(design$Z %*% result$fit$gamma) else rep(0, nrow(data))
    hazard <- numeric(nrow(data))
    truncated <- 0L
    for (i in seq_len(nrow(data))) {
        selected <- times > start[i] & if (left_limit) times < endpoint[i] else times <= endpoint[i]
        increment <- as.numeric(jumps[selected, , drop = FALSE] %*% design$X[i, ])
        # Additive hazard estimates need not be nonnegative. Use the positive
        # part of each fitted jump rather than allow survival above one.
        truncated <- truncated + sum(increment < -sqrt(.Machine$double.eps))
        hazard[i] <- sum(pmax(0, increment))
        if (result$model == "cox.aalen") {
            hazard[i] <- hazard[i] * exp(fixed[i])
        } else if (design$pz) {
            truncated <- truncated + as.integer(fixed[i] < -sqrt(.Machine$double.eps))
            hazard[i] <- hazard[i] + pmax(0, fixed[i]) * (endpoint[i] - start[i])
        }
    }
    if (truncated)
        warning("Negative fitted additive censoring-hazard increments were truncated at zero (", truncated,
                " prediction increments).", call. = FALSE)
    probabilities <- exp(-hazard)
    continuous_nuisance_predictions(list(predicted_values = probabilities), nrow(data), "Censoring", TRUE)
}
