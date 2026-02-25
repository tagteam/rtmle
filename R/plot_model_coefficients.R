### plot_nuisance.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 23 2026 (06:38) 
## Version: 
## Last-Updated: feb 23 2026 (10:41) 
##           By: Thomas Alexander Gerds
##     Update #: 10
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Plot extracted beta coefficients from nested intervention/outcome models
#'
#' @description
#' Extracts beta coefficients from a nested model object \code{x} (typically stored
#' under \code{x$models} with elements \code{time_0}, \dots, \code{time_K}) and
#' visualizes them as either (i) a dot plot with model/outcome identifiers on the
#' x-axis, or (ii) a single-panel "manhattan-style" plot in which all coefficients
#' from all selected models are shown in one graph.
#'
#' The function expects each relevant model entry to be a list containing at least
#' \code{$fit}, a sparse coefficient vector (a \code{Matrix} object of class
#' \code{dgCMatrix} with dimension \code{p x 1}) whose row names are coefficient
#' (term) names. Optional components \code{$formula} and \code{$warnings} are
#' ignored for plotting, but \code{$warnings} are returned in the output.
#'
#' @param x
#' An object containing nested fitted models. Must contain \code{x$models}, a
#' named list with elements like \code{time_0}, \code{time_1}, \dots. Each
#' \code{time_k} element is expected to be a list of model nodes (e.g.,
#' intervention protocol nodes, \code{censoring}, \code{outcome}), each of which
#' is itself a named list of model entries.
#'
#' @param node
#' Character scalar. For \code{plot_style = "by_outcome"}, the name of the node
#' within each \code{time_k} block from which to extract coefficients (e.g.
#' \code{"outcome"} or \code{"censoring"}).
#'
#' @param term
#' Optional filter for coefficient names. If \code{NULL}, all coefficient terms
#' are retained (subject to \code{include_intercept}). If a length-1 character
#' string, it is treated as a regular expression and matched using
#' \code{grepl()}. If a character vector of length > 1, it is treated as an
#' explicit set of term names to keep.
#'
#' @param include_intercept
#' Logical. If \code{FALSE} (default), removes \code{"(Intercept)"} from the
#' extracted coefficients when present.
#'
#' @param times
#' Optional time selection. If \code{NULL}, all \code{time_*} elements are used.
#' If numeric, values \code{0:K} are translated to names \code{time_0}, \dots.
#' If character, interpreted as names of \code{time_*} elements.
#'
#' @param outcomes
#' Optional filter for model/outcome names within a node. If \code{NULL}, all
#' are retained. If a length-1 character string, treated as a regular expression
#' matched using \code{grepl()}. If a character vector of length > 1, treated as
#' an explicit set of model names to keep.
#'
#' @param color_by
#' For \code{plot_style = "by_outcome"} only. Character scalar controlling point
#' coloring: \code{"time"} colors by time index; \code{"none"} uses a single
#' color.
#'
#' @param facet_by
#' For \code{plot_style = "by_outcome"} only. Character scalar controlling
#' faceting: \code{"term"} facets by coefficient term; \code{"time"} facets by
#' time index; \code{"none"} uses a single panel.
#'
#' @param point_alpha
#' Numeric in \code{[0,1]}. Point transparency passed to \code{geom_point()}.
#'
#' @param point_size
#' Numeric. Point size passed to \code{geom_point()}.
#'
#' @param plot_style
#' Character scalar. Selects plotting style:
#' \describe{
#'   \item{\code{"by_outcome"}}{Dot plot with model/outcome identifiers on the
#'   x-axis and beta coefficients on the y-axis (optionally colored/faceted).}
#'   \item{\code{"manhattan"}}{Single-panel plot showing all coefficients from
#'   all selected models. The x-axis corresponds to an ordering of models by time
#'   (\code{time_0}, \dots, \code{time_K}) and, within time, by node grouping as
#'   specified by \code{node_order}.}
#' }
#'
#' @param node_order
#' For \code{plot_style = "manhattan"} only. Character vector specifying the
#' within-time ordering of node groups. The special value \code{"protocol"}
#' expands to \code{protocol_nodes}. Typical order is
#' \code{c("protocol", "censoring", "outcome")}.
#'
#' @param protocol_nodes
#' For \code{plot_style = "manhattan"} only. Character vector naming the
#' intervention protocol nodes to include and order where \code{node_order}
#' contains \code{"protocol"} (e.g., \code{c("Placebo", "Lira")}).
#'
#' @param manhattan_color_by
#' For \code{plot_style = "manhattan"} only. Character scalar controlling point
#' coloring: \code{"none"} (single color), \code{"node_group"} (protocol vs
#' censoring vs outcome), \code{"time"} (time index), or \code{"term"}
#' (coefficient term).
#'
#' @param show_x_labels
#' For \code{plot_style = "manhattan"} only. Logical. If \code{TRUE}, shows the
#' full categorical labels for each model position on the x-axis. For large
#' numbers of models this may be slow and unreadable; by default labels are
#' hidden.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{data}}{A data frame in long format with one row per extracted
#'   coefficient. Contains at least \code{time}, \code{time_k}, \code{node},
#'   \code{node_group}, \code{outcome}, \code{term}, \code{beta}, and, for
#'   \code{plot_style = "manhattan"}, the ordering variables \code{x_cat} and
#'   \code{x_index}.}
#'   \item{\code{plot}}{A \code{ggplot} object.}
#'   \item{\code{warnings}}{A data frame of model warnings (if present), with
#'   columns \code{time}, \code{node}, \code{outcome}, and \code{warning}; or
#'   \code{NULL} if no warnings were found.}
#' }
#'
#' @details
#' Coefficients are extracted from \code{$fit} entries assumed to be sparse column
#' matrices (dimension \code{p x 1}). Term names are obtained from row names of
#' \code{$fit}. If no row names are present, generic names \code{beta_1}, \dots,
#' are assigned.
#'
#' In \code{plot_style = "manhattan"}, the x-axis represents an ordering of
#' distinct model identifiers defined by time, node grouping, node name and model
#' name. All coefficient terms for a given model are plotted at the same x
#' position, producing a vertical "stack" of points per model.
#'
#' @seealso
#' \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_point}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_bw labs theme element_text element_blank position_jitter facet_wrap
#' @importFrom Matrix as.matrix
#'
#' @export
plot_model_coefficients <- function(
  x,
  node = "outcome",
  term = NULL,                 # NULL = all terms; otherwise regex or character vector of term names
  include_intercept = FALSE,
  times = NULL,                # NULL = all; otherwise numeric like 0:K or names like "time_9"
  outcomes = NULL,             # NULL = all; otherwise regex or character vector of outcome/model names
  color_by = c("time", "none"),
  facet_by = c("term", "time", "none"),
  point_alpha = 0.8,
  point_size = 2.2,
  # NEW:
  plot_style = c("by_outcome", "manhattan"),
  # only used for manhattan:
  node_order = c("protocol", "censoring", "outcome"),
  protocol_nodes = c("Placebo", "Lira"),
  # manhattan aesthetics:
  manhattan_color_by = c("none", "node_group", "time", "term"),
  show_x_labels = FALSE
) {
  color_by <- match.arg(color_by)
  facet_by <- match.arg(facet_by)
  plot_style <- match.arg(plot_style)
  manhattan_color_by <- match.arg(manhattan_color_by)

  if (is.null(x$models) || !is.list(x$models)) stop("Expected x$models to be a list.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  # ---- helpers ----
  .is_sparse_colvec <- function(m) {
    inherits(m, "Matrix") && length(dim(m)) == 2 && dim(m)[2] == 1
  }
  .coef_from_fit <- function(fit) {
    if (!.is_sparse_colvec(fit)) return(NULL)
    v <- drop(Matrix::as.matrix(fit))
    nm <- rownames(fit)
    if (!is.null(nm) && length(nm) == length(v)) names(v) <- nm
    v
  }
  .matches <- function(x, pat_or_vec) {
    if (is.null(pat_or_vec)) return(rep(TRUE, length(x)))
    if (length(pat_or_vec) == 1L) grepl(pat_or_vec, x) else x %in% pat_or_vec
  }
  .normalize_time_selection <- function(time_names, times) {
    if (is.null(times)) return(time_names)
    if (is.numeric(times)) return(intersect(paste0("time_", times), time_names))
    intersect(as.character(times), time_names)
  }
  .parse_time_k <- function(tn) suppressWarnings(as.integer(sub("^time_", "", tn)))

  .nodes_to_traverse <- function(plot_style, node, node_order, protocol_nodes) {
    if (plot_style == "by_outcome") return(node)
    out <- character(0)
    for (z in node_order) {
      if (z == "protocol") out <- c(out, protocol_nodes) else out <- c(out, z)
    }
    unique(out)
  }
  .node_group <- function(node_name, protocol_nodes) {
    if (node_name %in% protocol_nodes) return("protocol")
    if (node_name == "censoring") return("censoring")
    if (node_name == "outcome") return("outcome")
    "other"
  }
  .extract_from_time_node <- function(time_block, time_name, node_name) {
    nb <- time_block[[node_name]]
    if (is.null(nb) || !is.list(nb)) return(NULL)

    model_names <- names(nb)
    if (is.null(model_names)) return(NULL)

    model_names <- model_names[.matches(model_names, outcomes)]

    rows <- list()
    warn_rows <- list()

    for (mn in model_names) {
      mobj <- nb[[mn]]
      if (!is.list(mobj) || is.null(mobj$fit)) next

      coefs <- .coef_from_fit(mobj$fit)
      if (is.null(coefs)) next

      if (!is.null(mobj$warnings)) {
        warn_rows[[length(warn_rows) + 1L]] <- data.frame(
          time = time_name,
          node = node_name,
          outcome = mn,
          warning = paste(mobj$warnings, collapse = " | "),
          stringsAsFactors = FALSE
        )
      }

      if (!include_intercept && !is.null(names(coefs))) {
        coefs <- coefs[names(coefs) != "(Intercept)"]
      }

      if (is.null(names(coefs))) names(coefs) <- paste0("beta_", seq_along(coefs))

      keep_term <- .matches(names(coefs), term)
      coefs <- coefs[keep_term]
      if (length(coefs) == 0) next

      rows[[length(rows) + 1L]] <- data.frame(
        time = time_name,
        node = node_name,
        outcome = mn,
        term = names(coefs),
        beta = as.numeric(coefs),
        stringsAsFactors = FALSE
      )
    }

    list(
      data = if (length(rows)) do.call(rbind, rows) else NULL,
      warnings = if (length(warn_rows)) do.call(rbind, warn_rows) else NULL
    )
  }

  # ---- choose time elements ----
  time_names <- names(x$models)
  if (is.null(time_names)) stop("x$models must be a named list (e.g., time_0, time_1, ...).")
  sel_times <- .normalize_time_selection(time_names, times)
  if (length(sel_times) == 0) stop("No matching time_* elements found in x$models for 'times'.")

  # ---- extract long table ----
  nodes_to_traverse <- .nodes_to_traverse(plot_style, node, node_order, protocol_nodes)

  all_rows <- list()
  all_warn <- list()

  for (tn in sel_times) {
    tb <- x$models[[tn]]
    if (!is.list(tb)) next
    for (nd in nodes_to_traverse) {
      rr <- .extract_from_time_node(tb, tn, nd)
      if (!is.null(rr$data)) all_rows[[length(all_rows) + 1L]] <- rr$data
      if (!is.null(rr$warnings)) all_warn[[length(all_warn) + 1L]] <- rr$warnings
    }
  }
  if (length(all_rows) == 0) stop("No coefficients found (check plot_style/node/times).")

  df <- do.call(rbind, all_rows)
  warnings_df <- if (length(all_warn)) do.call(rbind, all_warn) else NULL

  df$time_k <- .parse_time_k(df$time)
  if (all(is.na(df$time_k))) df$time_k <- df$time
  df$node_group <- vapply(df$node, .node_group, character(1), protocol_nodes = protocol_nodes)

  # -------------------------
  # PLOTTING
  # -------------------------
  if (plot_style == "by_outcome") {
    gg <- ggplot2::ggplot(df, ggplot2::aes(x = outcome, y = beta))

    if (color_by == "time") {
      gg <- gg + ggplot2::geom_point(
        ggplot2::aes(color = factor(time_k)),
        alpha = point_alpha, size = point_size,
        position = ggplot2::position_jitter(width = 0.15, height = 0)
      ) + ggplot2::labs(color = "time")
    } else {
      gg <- gg + ggplot2::geom_point(
        alpha = point_alpha, size = point_size,
        position = ggplot2::position_jitter(width = 0.15, height = 0)
      )
    }

    gg <- gg +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::labs(
        x = "Outcome (model name)",
        y = "Beta coefficient",
        title = paste0("Extracted coefficients (node = '", node, "')")
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.minor = ggplot2::element_blank()
      )

    if (facet_by == "term") gg <- gg + ggplot2::facet_wrap(~ term, scales = "free_y")
    if (facet_by == "time") gg <- gg + ggplot2::facet_wrap(~ time_k)

    return(invisible(list(data = df, plot = gg, warnings = warnings_df)))
  }

  # ---- manhattan style: SINGLE PANEL ----
  # You asked: one panel; x-axis = "outcome variables" ordered by time_0..time_K and within time by:
  # protocol variables, censoring, outcome. y-axis = all corresponding betas.
  #
  # We'll create a single categorical x-axis label per (time, node_group, node, outcome),
  # and then plot ALL betas (terms) vertically at that x position.

  # Order node_group by requested order
  df$node_group <- factor(df$node_group, levels = unique(c(node_order, "other")))

  # Build the x "outcome variable" label
  # (includes node so Placebo vs Lira stays distinct)
  df$x_cat <- paste0(
    "time_", df$time_k, " | ",
    as.character(df$node_group), " | ",
    df$node, " | ",
    df$outcome
  )

  # Determine ordering for x_cat
  ord <- with(df, order(
    ifelse(is.na(time_k), as.integer(factor(time)), time_k),
    node_group,
    node,
    outcome
  ))
  x_levels <- unique(df$x_cat[ord])
  df$x_cat <- factor(df$x_cat, levels = x_levels)

  # Manhattan-like numeric x index (better for huge axes)
  df$x_index <- as.integer(df$x_cat)

  # Aesthetics mapping
  aes_base <- ggplot2::aes(x = x_index, y = beta)
  if (manhattan_color_by == "time") {
    aes_base <- ggplot2::aes(x = x_index, y = beta, color = factor(time_k))
  } else if (manhattan_color_by == "node_group") {
    aes_base <- ggplot2::aes(x = x_index, y = beta, color = node_group)
  } else if (manhattan_color_by == "term") {
    aes_base <- ggplot2::aes(x = x_index, y = beta, color = term)
  }

  gg <- ggplot2::ggplot(df, aes_base) +
      ggplot2::geom_point(alpha = point_alpha, size = point_size,
                          ggplot2::aes(
                                       text = paste0(
                                           "Time: ", time_k,
                                           "<br>Node group: ", node_group,
                                           "<br>Node: ", node,
                                           "<br>Outcome: ", outcome,
                                           "<br>Term: ", term,
                                           "<br>Beta: ", signif(beta, 4)
                                       ))) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::labs(
                   x = "Models ordered by time and node group",
                   y = "Beta coefficient",
                   title = "All betas across all models (single-panel manhattan style)",
                   color = if (manhattan_color_by == "none") NULL else manhattan_color_by
               ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x  = if (show_x_labels)
                                      ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
                                  else
                                      ggplot2::element_blank()
               )

  if (show_x_labels) {
    # NOTE: for many x-levels this will be unreadable/slow.
    gg <- gg + ggplot2::scale_x_continuous(
      breaks = seq_along(x_levels),
      labels = x_levels
    )
  }

  invisible(list(data = df, plot = gg, warnings = warnings_df))
}

######################################################################
### plot_nuisance.R ends here
