#' Fit an rtmle analysis separately by stratum
#'
#' Selects participants by factor level, numeric value, or numeric interval
#' and refits the analysis in each group.
#'
#' @param x An \code{rtmle} object with prepared data, regimes, targets, and
#'   model formulas.
#' @param variable A baseline covariate supplied through
#'   \code{\link{add_baseline_data}} and present in \code{x$prepared_data},
#'   given as a symbol or a character string.
#' @param breaks An increasing vector of at least two interval boundaries for
#'   a numeric covariate. Use \code{-Inf} and \code{Inf} for unbounded interval boundaries.
#'   Intervals include their right endpoints and the first interval also includes
#'   its left endpoint. With \code{breaks = NULL}, factors and numeric covariates
#'   are grouped by their observed values.
#' @param ... Arguments passed to \code{\link{run_rtmle}}, such as
#'   \code{learner}, \code{estimator}, and \code{time_horizon}.
#' @details Empty factor levels and intervals are skipped. Participants with
#'   missing values or values outside the specified intervals are excluded with
#'   a warning. The function also warns when there are more than five nonempty
#'   strata or when a stratum contains fewer than 30 participants.
#'   The groups are used only to select participants.
#' @return An object of class \code{stratified_rtmle}, a named list containing
#'   one fitted \code{rtmle} object per stratum. The attributes \code{variable},
#'   \code{breaks} (when supplied), and \code{counts} give the column name,
#'   interval boundaries, and participant counts, respectively.
#' @seealso \code{\link{subset.rtmle}}, \code{\link{run_rtmle}},
#'   \code{\link{summary.stratified_rtmle}}, \code{\link{plot.stratified_rtmle}}
#' @examples
#' data(rtmle_object)
#'
#' # Fit separate analyses for the two values of the numeric indicator sex.
#' by_sex <- stratify.rtmle(rtmle_object, sex, time_horizon = 1:2)
#' print(by_sex)
#'
#' # Use break points to select age groups and retain age as a numeric covariate.
#' by_age <- stratify.rtmle(rtmle_object, "age", breaks = c(-Inf, 50, Inf),
#'                         time_horizon = 1:2)
#' summary(by_age)
#' @export
stratify <- function(x, ...) {
    UseMethod("stratify")
}

#' @rdname stratify
#' @export stratify.rtmle
#' @export
stratify.rtmle <- function(x, variable, breaks = NULL, ...) {
    if (!inherits(x, "rtmle")) {
        stop("x must be an rtmle object.", call. = FALSE)
    }
    data <- x$prepared_data
    id_name <- x$names$id
    if (NROW(data) == 0L || !id_name %in% names(data)) {
        stop("Prepare the data with prepare_rtmle_data() before stratifying.",
             call. = FALSE)
    }
    if (missing(variable)) {
        stop("Supply a variable in x$prepared_data.", call. = FALSE)
    }
    column <- substitute(variable)
    if (is.symbol(column) && as.character(column) %in% names(data)) {
        column <- as.character(column)
    } else {
        column <- tryCatch(variable, error = function(error) NULL)
    }
    if (!is.character(column) || length(column) != 1L || !column %in% names(data)) {
        stop("variable must be a column in x$prepared_data.", call. = FALSE)
    }
    if (!column %in% x$names$name_baseline_covariates) {
        stop("Stratifying variable ", sQuote(column), " must be a baseline covariate. ",
             "Supply baseline covariates with add_baseline_data() before preparing the data.",
             call. = FALSE)
    }

    values <- data[[column]]
    if (is.factor(values)) {
        if (!is.null(breaks)) {
            stop("Set breaks to NULL when stratifying by a factor.", call. = FALSE)
        }
        groups <- droplevels(values, exclude = NA)
    } else if (!is.numeric(values) || !is.null(dim(values))) {
        stop("The stratifying variable must be a factor or a numeric vector.",
             call. = FALSE)
    } else if (is.null(breaks)) {
        groups <- factor(values, levels = sort(unique(values)))
    } else {
        if (!is.numeric(breaks) || !is.null(dim(breaks)) || length(breaks) < 2L ||
            anyNA(breaks) || is.unsorted(breaks, strictly = TRUE)) {
            stop("For a numeric variable, supply at least two strictly increasing breaks.",
                 call. = FALSE)
        }
        groups <- cut(values, breaks = breaks, include.lowest = TRUE)
    }

    strata_ids <- split(data[[id_name]], groups, drop = TRUE)
    if (!length(strata_ids)) {
        stop("No participants fall into a stratum of ", sQuote(column), ".",
             call. = FALSE)
    }
    excluded <- sum(is.na(groups))
    if (excluded > 0L) {
        warning("Excluded ", excluded, " participant(s) with missing or out-of-range ",
                "values of ", sQuote(column), ".", call. = FALSE)
    }
    counts <- lengths(strata_ids)
    if (length(counts) > 5L) {
        warning("Stratifying by ", sQuote(column), " produced ", length(counts),
                " nonempty strata (more than 5).", call. = FALSE)
    }
    small_strata <- counts[counts < 30]
    if (length(small_strata) > 0L) {
        labels <- paste0(names(small_strata), " (n = ", small_strata, ")")
        warning("Small strata for ", sQuote(column), " (fewer than 30 participants): ", paste(labels, collapse = ", "),
                ".", call. = FALSE)
    }
    fits <- lapply(seq_along(strata_ids), function(i) {
        tryCatch(
            run_rtmle(subset(x, ids = strata_ids[[i]]), ...),
            error = function(error) {
                stop("Could not fit stratum ", sQuote(names(strata_ids)[[i]]),
                     " of ", sQuote(column), ": ", conditionMessage(error),
                     call. = FALSE)
            }
        )
    })
    names(fits) <- names(strata_ids)
    structure(fits, class = "stratified_rtmle", variable = column,
              breaks = breaks, counts = counts)
}
