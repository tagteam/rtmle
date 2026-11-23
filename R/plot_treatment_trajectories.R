### plot_treatment_trajectories.R ---
#----------------------------------------------------------------------
## Plot individual treatment trajectories in continuous and discrete time.
#----------------------------------------------------------------------
#' Plot individual treatment trajectories
#'
#' Shows subject-level treatment histories from an \code{rtmle} object. When
#' long-format treatment data are available, the continuous-time panel uses
#' the histories in \code{x$long_data$timevar_data}. Point times use the
#' configured time column or \code{date}. When discretized
#' data are available, the discrete-time panel uses \code{x$data$timevar_data}
#' or \code{x$prepared_data}. If both representations are present, the panels
#' are displayed side by side.
#'
#' For multiple treatment variables, colour identifies the joint treatment
#' state (for example, \code{A=1, B=0}), preserving one horizontal trajectory
#' per subject. In the continuous-time panel, observed point-format values are
#' overlaid as points; exposure-interval starts and ends are marked as points.
#' Legend entries are sorted, with states containing \code{Missing} last.
#'
#' @param x An object of class \code{"rtmle"} containing long-format and/or
#'   wide-format treatment data.
#' @param treatment_variables Character vector naming the treatment variables.
#'   By default, names in \code{x$names$treatment_options} (or the treatment
#'   variables of a stored regime) are used. Supplying this argument is useful
#'   before a regime has been added to an object.
#' @return A \code{ggplot2} object. The y-axis contains subjects numbered from
#'   1 to the number of subjects, but subject labels are intentionally hidden.
#' @seealso \code{\link{add_long_data}}, \code{\link{discretize_data}},
#'   \code{\link{add_wide_data}}, \code{\link{prepare_rtmle_data}}
#' @examples
#' data(rtmle_object)
#' p <- plot_treatment_trajectories(rtmle_object)
#' p
#' @importFrom ggplot2 aes facet_grid geom_point geom_segment geom_vline ggplot labs
#'   scale_colour_manual scale_x_continuous scale_y_continuous theme theme_bw
#' @export
plot_treatment_trajectories <- function(x, treatment_variables = NULL) {

    start <- end <- subject <- state <- xintercept <- time <- NULL


    if (!inherits(x, "rtmle")) {
        stop("Argument x must be an object of class 'rtmle'.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' must be installed to use plot_treatment_trajectories().")
    }

    id <- x$names$id
    if (is.null(id) || length(id) != 1L || !nzchar(id)) {
        stop("The rtmle object does not contain a valid subject identifier name.")
    }

    helper_plot_long_table <- function(variable) {
        if (is.null(x$long_data) || is.null(x$long_data$timevar_data)) {
            return(NULL)
        }
        timevar_data <- x$long_data$timevar_data
        if (inherits(timevar_data, "data.frame")) {
            return(timevar_data)
        }
        if (is.list(timevar_data) && variable %in% names(timevar_data)) {
            return(timevar_data[[variable]])
        }
        NULL
    }

    helper_plot_wide_table <- function(variable) {
        # Continuous-time prepared columns are indexed by event order.
        if (isTRUE(x$continuous_time)) return(NULL)
        candidates <- list()
        timevar_data <- if (!is.null(x$data)) x$data$timevar_data else NULL

        if (inherits(timevar_data, "data.frame")) {
            candidates[[length(candidates) + 1L]] <- timevar_data
        } else if (is.list(timevar_data) && variable %in% names(timevar_data)) {
            candidates[[length(candidates) + 1L]] <- timevar_data[[variable]]
        }
        if (!is.null(x$prepared_data)) {
            candidates[[length(candidates) + 1L]] <- x$prepared_data
        }

        prefix <- paste0(
            "^", gsub("([.|()\\^{}+$*?]|\\\\)", "\\\\\\1", variable),
            "_[0-9]+$"
        )
        for (candidate in candidates) {
            if (!inherits(candidate, "data.frame") || !(id %in% names(candidate))) {
                next
            }
            columns <- grep(prefix, names(candidate), value = TRUE)
            if (length(columns) > 0L) {
                suffix <- sub("^.*_", "", columns)
                return(list(
                    data = candidate,
                    columns = columns,
                    nodes = stats::setNames(columns, suffix)
                ))
            }
        }
        NULL
    }

    helper_plot_time_column <- function(data, fallback = "date") {
        columns <- intersect(c(x$names$time, fallback), setdiff(names(data), id))
        if (length(columns) > 0L) columns[[1L]] else NULL
    }

    helper_plot_time_values <- function(values, ids, label) {
        if (inherits(values, "Date") || inherits(values, "POSIXt")) {
            start_name <- x$names$start_followup_date
            baseline <- if (!is.null(x$data)) x$data$baseline_data else NULL
            if (is.null(start_name) || length(start_name) == 0L ||
                !inherits(baseline, "data.frame") ||
                !(id %in% names(baseline)) || !(start_name %in% names(baseline))) {
                stop("Cannot plot ", label,
                     " because its dates are calendar dates and the object does not contain ",
                     "a subject-specific start_followup_date column.")
            }
            starts <- baseline[[start_name]][
                match(as.character(ids), as.character(baseline[[id]]))
            ]
            if (anyNA(starts)) {
                stop("The start_followup_date values needed to plot ", label,
                     " contain missing values.")
            }
            return(as.numeric(values - starts))
        }
        as.numeric(values)
    }

    helper_plot_value_label <- function(value) {
        if (length(value) == 0L || is.na(value[[1L]])) {
            return("Missing")
        }
        value <- value[[1L]]
        if (is.logical(value)) {
            return(if (isTRUE(value)) "1" else "0")
        }
        as.character(value)
    }

    helper_plot_state <- function(values, variables) {
        paste(
            paste0(
                variables, "=",
                vapply(values, helper_plot_value_label, character(1L))
            ),
            collapse = ", "
        )
    }

    helper_plot_long_values <- function(table_info, subject_rows, time) {
        lapply(names(table_info), function(variable) {
            info <- table_info[[variable]]
            rows <- subject_rows[[variable]]
            if (length(rows) == 0L) {
                return(if (info$interval) 0 else NA)
            }
            if (info$interval) {
                starts <- info$starts[rows]
                ends <- info$ends[rows]
                exposed <- any(
                    is.finite(starts) & starts <= time &
                        (is.na(ends) | !is.finite(ends) | ends > time)
                )
                return(if (exposed) 1 else 0)
            }
            times <- info$times[rows]
            observed <- which(
                !is.na(times) & is.finite(times) & times <= time &
                    !is.na(info$values[rows])
            )
            if (length(observed) == 0L) {
                return(NA)
            }
            last_observation <- observed[[which.max(times[observed])]]
            info$values[rows][last_observation]
        })
    }

    helper_plot_long_segments <- function(variables, subject_ids, subject_rank,
                                          end_followup) {
        table_info <- list()
        for (variable in variables) {
            data <- helper_plot_long_table(variable)
            if (!inherits(data, "data.frame") || NROW(data) == 0L ||
                !(id %in% names(data))) {
                next
            }
            interval_format <- all(c("start_date", "end_date") %in% names(data))
            time_column <- helper_plot_time_column(data)
            point_format <- !is.null(time_column) && "value" %in% names(data)
            if (!interval_format && !point_format) {
                stop("Long-format treatment data for '", variable,
                     "' must contain an id with exposure intervals or point times and values.")
            }
            data_ids <- as.character(data[[id]])
            if (interval_format) {
                starts <- helper_plot_time_values(
                    data[["start_date"]], data_ids,
                    paste0("long_data$timevar_data$", variable)
                )
                ends <- helper_plot_time_values(
                    data[["end_date"]], data_ids,
                    paste0("long_data$timevar_data$", variable)
                )
                times <- NULL
                values <- NULL
            } else {
                times <- helper_plot_time_values(
                    data[[time_column]], data_ids,
                    paste0("long_data$timevar_data$", variable)
                )
                values <- data[["value"]]
                starts <- NULL
                ends <- NULL
            }
            table_info[[variable]] <- list(
                ids = data_ids,
                interval = interval_format,
                starts = starts,
                ends = ends,
                times = times,
                values = values
            )
        }
        if (length(table_info) == 0L) {
            return(NULL)
        }

        out <- vector("list", 0L)
        points <- vector("list", 0L)
        for (subject_index in seq_along(subject_ids)) {
            subject_id <- subject_ids[[subject_index]]
            end <- end_followup[[subject_index]]
            if (is.na(end) || end <= 0) {
                next
            }
            breakpoints <- c(0, end)
            subject_rows <- list()
            for (variable in names(table_info)) {
                info <- table_info[[variable]]
                rows <- which(info$ids == subject_id)
                subject_rows[[variable]] <- rows
                if (length(rows) == 0L) {
                    next
                }
                if (info$interval) {
                    starts <- info$starts[rows]
                    ends <- info$ends[rows]
                    breakpoints <- c(
                        breakpoints,
                        starts[is.finite(starts) & starts > 0 & starts < end],
                        ends[is.finite(ends) & ends > 0 & ends < end]
                    )
                } else {
                    times <- info$times[rows]
                    breakpoints <- c(
                        breakpoints,
                        times[is.finite(times) & times > 0 & times < end]
                    )
                }
            }
            breakpoints <- sort(unique(breakpoints[is.finite(breakpoints)]))
            if (length(breakpoints) < 2L) {
                next
            }

            for (segment_index in seq_len(length(breakpoints) - 1L)) {
                start <- breakpoints[[segment_index]]
                finish <- breakpoints[[segment_index + 1L]]
                if (finish <= start) {
                    next
                }
                midpoint <- (start + finish) / 2
                values <- helper_plot_long_values(
                    table_info, subject_rows, midpoint
                )
                state <- helper_plot_state(values, names(table_info))
                out[[length(out) + 1L]] <- data.frame(
                    panel = "Continuous time",
                    start = start,
                    end = finish,
                    subject = subject_rank[[subject_index]],
                    state = state,
                    stringsAsFactors = FALSE
                )
            }

            ## Add the raw observed values as points.  Point-format histories
            ## contribute one point per observed value.  For interval-format
            ## histories the start and end of each observed exposure interval
            ## are marked, including the clipped baseline/follow-up boundary.
            for (variable in names(table_info)) {
                info <- table_info[[variable]]
                rows <- subject_rows[[variable]]
                if (length(rows) == 0L) {
                    next
                }
                if (info$interval) {
                    starts <- info$starts[rows]
                    ends <- info$ends[rows]
                    for (row_index in seq_along(rows)) {
                        interval_start <- starts[[row_index]]
                        interval_end <- ends[[row_index]]
                        if (is.finite(interval_start) &&
                            (is.na(interval_end) ||
                             !is.finite(interval_end) ||
                             interval_end > 0)) {
                            point_time <- max(0, interval_start)
                            if (point_time <= end &&
                                (is.na(interval_end) ||
                                 !is.finite(interval_end) ||
                                 interval_end > point_time)) {
                                values <- helper_plot_long_values(
                                    table_info, subject_rows, point_time
                                )
                                points[[length(points) + 1L]] <- data.frame(
                                    panel = "Continuous time",
                                    time = point_time,
                                    subject = subject_rank[[subject_index]],
                                    state = helper_plot_state(
                                        values, names(table_info)
                                    ),
                                    stringsAsFactors = FALSE
                                )
                            }
                        }
                        if (is.finite(interval_end) &&
                            interval_end >= 0 && interval_end <= end &&
                            (is.na(interval_start) ||
                             interval_end >= interval_start)) {
                            values <- helper_plot_long_values(
                                table_info, subject_rows, interval_end
                            )
                            points[[length(points) + 1L]] <- data.frame(
                                panel = "Continuous time",
                                time = interval_end,
                                subject = subject_rank[[subject_index]],
                                state = helper_plot_state(
                                    values, names(table_info)
                                ),
                                stringsAsFactors = FALSE
                            )
                        }
                    }
                } else {
                    times <- info$times[rows]
                    observed <- which(
                        !is.na(times) & is.finite(times) &
                            times >= 0 & times <= end &
                            !is.na(info$values[rows])
                    )
                    for (row_index in observed) {
                        point_time <- times[[row_index]]
                        values <- helper_plot_long_values(
                            table_info, subject_rows, point_time
                        )
                        points[[length(points) + 1L]] <- data.frame(
                            panel = "Continuous time",
                            time = point_time,
                            subject = subject_rank[[subject_index]],
                            state = helper_plot_state(
                                values, names(table_info)
                            ),
                            stringsAsFactors = FALSE
                        )
                    }
                }
            }
        }
        if (length(out) == 0L) {
            return(NULL)
        }
        list(
            segments = do.call(rbind, out),
            points = if (length(points) > 0L) {
                do.call(rbind, points)
            } else {
                NULL
            }
        )
    }

    helper_plot_discrete_segments <- function(variables, subject_ids, subject_rank,
                                              end_followup, grid) {
        table_info <- lapply(variables, helper_plot_wide_table)
        names(table_info) <- variables
        table_info <- table_info[!vapply(table_info, is.null, logical(1L))]
        if (length(table_info) == 0L) {
            return(NULL)
        }

        node_numbers <- sort(unique(unlist(lapply(table_info, function(info) {
            as.integer(names(info$nodes))
        }), use.names = FALSE)))
        node_numbers <- node_numbers[
            is.finite(node_numbers) & node_numbers >= 0 & node_numbers <= length(grid) - 2L
        ]
        if (length(node_numbers) == 0L) {
            return(NULL)
        }

        out <- vector("list", 0L)
        for (subject_index in seq_along(subject_ids)) {
            subject_id <- subject_ids[[subject_index]]
            end <- end_followup[[subject_index]]
            if (is.na(end) || end <= 0) {
                next
            }
            for (node in node_numbers) {
                start <- grid[[node + 1L]]
                finish <- min(grid[[node + 2L]], end)
                if (finish <= start) {
                    next
                }
                values <- lapply(names(table_info), function(variable) {
                    info <- table_info[[variable]]
                    column <- info$nodes[[as.character(node)]]
                    if (is.null(column)) {
                        return(NA)
                    }
                    row <- match(subject_id, as.character(info$data[[id]]))
                    if (is.na(row)) {
                        return(NA)
                    }
                    info$data[[column]][[row]]
                })
                state <- paste(
                    paste0(
                        names(table_info), "=",
                        vapply(values, function(value) {
                            if (length(value) == 0L || is.na(value[[1L]])) {
                                return("Missing")
                            }
                            value <- value[[1L]]
                            if (is.logical(value)) {
                                return(if (isTRUE(value)) "1" else "0")
                            }
                            as.character(value)
                        }, character(1L))
                    ),
                    collapse = ", "
                )
                out[[length(out) + 1L]] <- data.frame(
                    panel = "Discrete time",
                    start = start,
                    end = finish,
                    subject = subject_rank[[subject_index]],
                    state = state,
                    stringsAsFactors = FALSE
                )
            }
        }
        if (length(out) == 0L) NULL else do.call(rbind, out)
    }

    if (!is.null(treatment_variables)) {
        if (!is.character(treatment_variables) ||
            length(treatment_variables) == 0L ||
            anyNA(treatment_variables) ||
            any(!nzchar(treatment_variables))) {
            stop("Argument treatment_variables must be a non-empty character vector.")
        }
        variables <- unique(treatment_variables)
    } else {
        variables <- character()
        if (!is.null(x$names$treatment_options)) {
            variables <- c(variables, names(x$names$treatment_options))
        }
        if (!is.null(x$regimes) && is.list(x$regimes)) {
            variables <- c(
                variables,
                unlist(lapply(x$regimes, function(regime) regime$treatment_variables),
                       use.names = FALSE)
            )
        }
        if (!is.null(x$names$treatment_variables)) {
            variables <- c(variables, x$names$treatment_variables)
        }
        variables <- unique(variables[!is.na(variables) & nzchar(variables)])

        if (length(variables) == 0L && isTRUE(x$continuous_time)) {
            variables <- "A"
        }
        if (length(variables) == 0L) {
            long_names <- if (!is.null(x$long_data$timevar_data) &&
                              is.list(x$long_data$timevar_data) &&
                              !inherits(x$long_data$timevar_data, "data.frame")) {
                names(x$long_data$timevar_data)
            } else {
                character()
            }
            long_names <- long_names[vapply(
                long_names,
                function(variable) {
                    data <- helper_plot_long_table(variable)
                    if (!inherits(data, "data.frame")) return(FALSE)
                    if (all(c("start_date", "end_date") %in% names(data))) {
                        return(TRUE)
                    }
                    if (is.null(helper_plot_time_column(data)) || !"value" %in% names(data)) {
                        return(FALSE)
                    }
                    values <- data[["value"]]
                    length(unique(as.character(values[!is.na(values)]))) <= 2L
                },
                logical(1L)
            )]
            wide_data <- if (!is.null(x$data)) x$data$timevar_data else NULL
            wide_tables <- if (inherits(wide_data, "data.frame")) {
                list(wide_data)
            } else if (is.list(wide_data)) {
                wide_data
            } else {
                list()
            }
            wide_tables <- c(wide_tables, list(x$prepared_data))
            wide_names <- unlist(lapply(wide_tables, function(table) {
                if (!inherits(table, "data.frame")) return(character())
                sub("_[0-9]+$", "", grep("_[0-9]+$", names(table), value = TRUE))
            }), use.names = FALSE)
            variables <- unique(intersect(long_names, wide_names))
            if (length(variables) == 0L && length(long_names) > 0L) {
                variables <- long_names
            }
            if (length(variables) == 0L && length(unique(wide_names)) == 1L) {
                variables <- unique(wide_names)
            }
        }
        if (length(variables) == 0L) {
            stop("Cannot determine treatment variables. Supply them with treatment_variables = c(...).")
        }
    }
    grid <- x$time_grid_scale
    if (is.null(grid) || length(grid) == 0L) grid <- x$time_grid
    if (is.null(grid) || length(grid) == 0L) grid <- c(0, 1)
    grid <- as.numeric(grid)
    if (length(grid) < 2L || any(!is.finite(grid)) || any(diff(grid) <= 0)) {
        stop("x$time_grid_scale must contain at least two increasing finite time points.")
    }

    id_tables <- list()
    if (!is.null(x$data)) {
        id_tables <- c(id_tables, list(
            x$data$baseline_data,
            x$data$outcome_data,
            x$data$censored_data,
            x$data$competing_data
        ))
        if (inherits(x$data$timevar_data, "data.frame")) {
            id_tables <- c(id_tables, list(x$data$timevar_data))
        } else if (is.list(x$data$timevar_data)) {
            id_tables <- c(id_tables, x$data$timevar_data)
        }
    }
    if (!is.null(x$long_data)) {
        id_tables <- c(id_tables, list(
            x$long_data$outcome_data,
            x$long_data$censored_data,
            x$long_data$competing_data
        ))
        if (inherits(x$long_data$timevar_data, "data.frame")) {
            id_tables <- c(id_tables, list(x$long_data$timevar_data))
        } else if (is.list(x$long_data$timevar_data)) {
            id_tables <- c(id_tables, x$long_data$timevar_data)
        }
    }
    id_tables <- c(id_tables, list(x$prepared_data, x$followup))
    subject_ids <- unique(unlist(lapply(id_tables, function(data) {
        if (!inherits(data, "data.frame") || !(id %in% names(data))) {
            return(character())
        }
        ids <- as.character(data[[id]])
        ids[!is.na(ids)]
    }), use.names = FALSE))
    if (length(subject_ids) == 0L) {
        stop("Cannot plot treatment trajectories because the object contains no subject ids.")
    }
    subject_rank <- seq_along(subject_ids)
    names(subject_rank) <- subject_ids

    end_followup <- rep(NA_real_, length(subject_ids))
    names(end_followup) <- subject_ids
    long_data <- x$long_data
    event_names <- intersect(
        c("competing_data", "outcome_data", "censored_data"),
        if (is.null(long_data)) character() else names(long_data)
    )
    for (event_name in event_names) {
        event_data <- long_data[[event_name]]
        if (!inherits(event_data, "data.frame") || NROW(event_data) == 0L ||
            !(id %in% names(event_data))) {
            next
        }
        time_column <- helper_plot_time_column(event_data, fallback = c("date", "time"))
        if (is.null(time_column)) next

        event_ids <- as.character(event_data[[id]])
        event_times <- helper_plot_time_values(
            event_data[[time_column]], event_ids,
            paste0("long_data$", event_name)
        )
        valid <- !is.na(event_times) & is.finite(event_times) & !is.na(event_ids)
        if (!any(valid)) next

        event_min <- tapply(event_times[valid], event_ids[valid], min)
        positions <- match(names(event_min), subject_ids)
        keep <- !is.na(positions)
        positions <- positions[keep]
        event_min <- as.numeric(event_min[keep])
        replace <- is.na(end_followup[positions]) | event_min < end_followup[positions]
        end_followup[positions[replace]] <- event_min[replace]
    }

    ## If exact event dates are not available, use the discrete follow-up
    ## information when it is present. last_interval counts the interval
    ## containing the last observed at-risk time, hence the +2 index here.
    followup <- x$followup
    if (anyNA(end_followup) && inherits(followup, "data.frame") &&
        all(c(id, "last_interval") %in% names(followup))) {
        followup_ids <- as.character(followup[[id]])
        interval <- as.numeric(followup[["last_interval"]])
        positions <- match(subject_ids, followup_ids)
        valid <- !is.na(positions) & !is.na(interval[positions])
        if (any(valid)) {
            grid_position <- pmin(
                pmax(floor(interval[positions[valid]]) + 2, 2),
                length(grid)
            )
            end_followup[valid] <- grid[grid_position]
        }
    }
    end_followup[is.na(end_followup) | !is.finite(end_followup)] <- max(grid)
    end_followup <- pmax(end_followup, 0)

    long_variables <- variables[vapply(
        variables,
        function(variable) {
            data <- helper_plot_long_table(variable)
            inherits(data, "data.frame") && NROW(data) > 0L
        },
        logical(1L)
    )]
    wide_variables <- variables[vapply(
        variables,
        function(variable) !is.null(helper_plot_wide_table(variable)),
        logical(1L)
    )]
    if (length(long_variables) == 0L && length(wide_variables) == 0L) {
        stop("No treatment data were found for the requested treatment_variables.")
    }

    continuous_result <- if (length(long_variables) > 0L) {
        helper_plot_long_segments(
            long_variables, subject_ids, subject_rank, end_followup
        )
    } else {
        NULL
    }
    continuous <- if (is.null(continuous_result)) {
        NULL
    } else {
        continuous_result$segments
    }
    points <- if (is.null(continuous_result)) {
        NULL
    } else {
        continuous_result$points
    }
    discrete <- if (length(wide_variables) > 0L) {
        helper_plot_discrete_segments(
            wide_variables, subject_ids, subject_rank, end_followup, grid
        )
    } else {
        NULL
    }
    if (is.null(continuous) && is.null(discrete)) {
        stop("The treatment data do not contain any plot-able trajectories.")
    }

    segments <- do.call(rbind, Filter(NROW, list(continuous, discrete)))
    panel_levels <- unique(as.character(segments$panel))
    segments$panel <- factor(segments$panel, levels = panel_levels)
    if (!is.null(points)) {
        points$panel <- factor(points$panel, levels = panel_levels)
    }
    state_values <- c(
        as.character(segments$state),
        if (is.null(points)) character() else as.character(points$state)
    )
    all_states <- unique(state_values)
    missing_state <- all_states == "Missing" |
        grepl("=Missing($|, )", all_states)
    states <- c(
        sort(all_states[!missing_state]),
        sort(all_states[missing_state])
    )
    if (length(states) == 0L) {
        stop("The treatment data do not contain any plot-able states.")
    }
    segments$state <- factor(segments$state, levels = states)
    if (!is.null(points)) {
        points$state <- factor(points$state, levels = states)
    }
    nonmissing <- states[!grepl("=Missing($|, )", states) & states != "Missing"]
    palette <- if (length(nonmissing) > 0L) {
        if (length(nonmissing) <= 8L) {
            c("#000000", "#E69F00", "#56B4E9", "#009E73",
              "#D55E00", "#0072B2", "#CC79A7", "#F0E442")[seq_along(nonmissing)]
        } else if (exists("hcl.colors", envir = asNamespace("grDevices"),
                          inherits = FALSE)) {
            grDevices::hcl.colors(length(nonmissing), palette = "Dark 3")
        } else {
            grDevices::rainbow(length(nonmissing), start = 0, end = 0.85)
        }
    } else {
        character()
    }
    names(palette) <- nonmissing
    missing_states <- setdiff(states, nonmissing)
    if (length(missing_states) > 0L) {
        missing_palette <- rep("#BDBDBD", length(missing_states))
        names(missing_palette) <- missing_states
        palette <- c(palette, missing_palette)
    }
    palette <- palette[states]

    vertical <- if ("Discrete time" %in% panel_levels) {
        data.frame(
            panel = factor("Discrete time", levels = panel_levels),
            xintercept = grid,
            stringsAsFactors = FALSE
        )
    } else {
        NULL
    }

    time_unit <- if (is.null(x$time_unit)) "Time" else x$time_unit
    p <- ggplot2::ggplot(segments) +
        ggplot2::geom_segment(
            ggplot2::aes(x = start, xend = end, y = subject, yend = subject,
                         colour = state),
            linewidth = 0.7,
            lineend = "butt"
        )
    if (!is.null(vertical)) {
        p <- ggplot2::ggplot(segments) +
            ggplot2::geom_vline(
                data = vertical,
                ggplot2::aes(xintercept = xintercept),
                inherit.aes = FALSE,
                colour = "grey75",
                linewidth = 0.35
            ) +
            ggplot2::geom_segment(
                ggplot2::aes(x = start, xend = end, y = subject, yend = subject,
                             colour = state),
                linewidth = 0.7,
                lineend = "butt"
            )
    }
    if (!is.null(points)) {
        p <- p + ggplot2::geom_point(
            data = points,
            mapping = ggplot2::aes(x = time, y = subject, colour = state),
            inherit.aes = FALSE,
            size = 2.1
        )
    }
    p <- p +
        ggplot2::scale_colour_manual(
            values = palette,
            drop = FALSE,
            name = "Exposure"
        ) +
        ggplot2::scale_x_continuous(
            breaks = grid,
            labels = if (length(x$time_grid_labels) == length(grid)) {
                x$time_grid_labels
            } else {
                grid
            },
            expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(1, length(subject_ids) + 1L),
            breaks = NULL,
            labels = NULL,
            expand = c(0, 0)
        ) +
        ggplot2::labs(x = time_unit, y = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )
    ## Keep the requested y scale (1 to N + 1), but allow the marker at the
    ## first subject position to extend into the lower plot margin.
    p <- p + ggplot2::coord_cartesian(clip = "off")
    if (length(panel_levels) > 1L) {
        p <- p + ggplot2::facet_grid(. ~ panel, scales = "free_x", space = "free_x")
    }
    p
}

######################################################################
### plot_treatment_trajectories.R ends here
