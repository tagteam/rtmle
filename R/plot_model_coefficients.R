### plot_nuisance.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 23 2026 (06:38) 
## Version: 
## Last-Updated: mar 26 2026 (15:49) 
##           By: Thomas Alexander Gerds
##     Update #: 133
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Plot beta coefficients from intervention/censoring/outcome models
#'
#' @description
#' Extracts beta coefficients from nuisance parameter models (stored
#' under \code{x$models} with elements \code{time_0}, \dots, \code{time_K}) and
#' visualizes them as either (i) a dot plot with model/outcome identifiers on the
#' x-axis, or (ii) a single-panel "manhattan-style" plot in which all coefficients
#' from all selected models are shown in one graph.
#'
#' @param x
#' An object fitted with \link{run_rtmle} containing fitted nuisance parameter models. 
#'
#' @param nodes
#' Character vector. For \code{plot_style = "by_outcome"}, the name of the node
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
#' @param time_horizon
#' Optional selection of time_horizon. If \code{NULL}, \code{max(x$time_grid)} is used.
#'
#' @param protocol
#' Optional selection of protocol. If \code{NULL}, \code{names(x$protocols)[[1]])} is used. 
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
#'   (\code{time_0}, \dots, \code{time_K}) and, within time.
#' }}
#'
#' @param manhattan_color_by
#' For \code{plot_style = "manhattan"} only. Character scalar controlling point
#' coloring: \code{"none"} (single color) (protocol vs
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
#'   \code{outcome}, \code{term}, \code{beta}, and, for
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
#' @examples
#' beta <- matrix(c(-1, 0.4), ncol = 1,
#'                dimnames = list(c("(Intercept)", "A_0"), "Estimate"))
#' x <- list(
#'   learner = list(name = "glm", fun = list("learn_glm")),
#'   times = 0:2,
#'   run_time_horizons = 2,
#'   time_grid = 0:2,
#'   intervention_nodes = 0:1,
#'   protocols = list(Always_A = list()),
#'   names = list(outcome = "Y"),
#'   models = list(
#'     time_0 = list(outcome = list(
#'       Y_1 = list(fit = list(Always_A = list(sequence_time_2 = beta))))),
#'     time_1 = list(outcome = list(
#'       Y_2 = list(fit = list(Always_A = list(sequence_time_2 = beta)))))))
#' class(x) <- "rtmle"
#' out <- plot_model_coefficients(x, time_horizon = 2, protocol = "Always_A",
#'                                nodes = "outcome", plot_style = "by_outcome")
#' names(out)
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_bw labs theme element_text element_blank position_jitter facet_wrap
#' @importFrom Matrix as.matrix
#' 
#' @export
plot_model_coefficients <- function(x,
                                    time_horizon,
                                    protocol,
                                    times = NULL,
                                    nodes,
                                    term = NULL,
                                    include_intercept = FALSE,
                                    outcomes = NULL,
                                    color_by = c("time", "none"),
                                    facet_by = c("term", "time", "none"),
                                    point_alpha = 0.8,
                                    point_size = 2.2,
                                    plot_style = c("manhattan", "by_outcome"),
                                    manhattan_color_by = c("node", "time", "term", "none"),
                                    show_x_labels = FALSE) {
    if (!(x$learner$fun[[1]] %chin% c("learn_glmnet","learn_glm"))){
        stop("Can only plot regression coefficients when fitter is either learn_glm or learn_glmnet")
    }
    outcome <- time_k <- x_index <- terms <- node <- x_cat <- NULL
    color_by <- match.arg(color_by)
    facet_by <- match.arg(facet_by)
    plot_style <- match.arg(plot_style)
    manhattan_color_by <- match.arg(manhattan_color_by)

    if (missing(time_horizon)){
        time_horizon <- max(x$run_time_horizons)
    }else{
        stopifnot(length(time_horizon) == 1)
        stopifnot(time_horizon %in% x$time_grid)
    }

    if (missing(protocol)){
        protocol <- names(x$protocols)[[1]]
    }else{
        stopifnot(length(protocol) == 1)
        stopifnot(protocol %in% names(x$protocols))
    }


    if (is.null(x$models) || !is.list(x$models)) stop("Expected x$models to be a list.")
    if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
    
    # ---- choose time elements ----
    time_names <- names(x$models)
    if (is.numeric(times)){
        times <- paste0("time_",times)
        selected_times <- intersect(time_names, times)
    }else{
        selected_times <- paste0("time_",x$time_grid[x$time_grid<time_horizon])
    }
    if (length(selected_times) == 0) stop("No matching time_* elements found in x$models for 'times'.")
    
    # ---- extract long table ----
    df <- coef(x)
    df[,time_k := suppressWarnings(as.integer(sub("^time_", "", time)))]
    if (!include_intercept){
        df <- df[terms != "(Intercept)"]
    }
    
    if (length(selected_times)>0){
        df <- df[time %in% selected_times]
    }


    if (!is.null(term)){
        df <- df[intersect(term,terms)]
    }
    ## if (plot_style == "by_outcome") {
        if (!missing(nodes)){
            stopifnot(all(nodes %in% df$node))
            df <- df[node %in% nodes]
        }
    ## }

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
                         title = ""
                     ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                         axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                         panel.grid.minor = ggplot2::element_blank()
                     )

        if (facet_by == "term") gg <- gg + ggplot2::facet_wrap(~ term, scales = "free_y")
        if (facet_by == "time") gg <- gg + ggplot2::facet_wrap(~ time_k)

        return(invisible(list(data = df, plot = gg)))
    }

    # ---- manhattan style: SINGLE PANEL ----
    # You asked: one panel; x-axis = "outcome variables" ordered by time_0..time_K and within time by:
    # protocol variables, censoring, outcome. y-axis = all corresponding betas.
    #
    # We'll create a single categorical x-axis label per (time, node, outcome),
    # and then plot ALL betas (terms) vertically at that x position.

    # Build the x "outcome variable" label
    # (includes node so Placebo vs Lira stays distinct)
    df[,x_cat := paste0("time_",time_k," | ",node," | ",outcome)]

    # Determine ordering for x_cat
    setkey(df,time,node,outcome)
    x_levels <- unique(df$x_cat)
    df[,x_cat := factor(x_cat, levels = x_levels)]
    # Manhattan-like numeric x index (better for huge axes)
    df[,x_index := as.integer(x_cat)]

    # Aesthetics mapping
    aes_base <- ggplot2::aes(x = x_index, y = beta)
    if (manhattan_color_by == "time") {
        aes_base <- ggplot2::aes(x = x_index, y = beta, color = factor(time_k))
    } else if (manhattan_color_by == "node") {
        aes_base <- ggplot2::aes(x = x_index, y = beta, color = node)
    } else if (manhattan_color_by == "terms") {
        aes_base <- ggplot2::aes(x = x_index, y = beta, color = terms)
    }

    gg <- ggplot2::ggplot(df, aes_base)
    suppressWarnings(gg <- gg+ ggplot2::geom_point(alpha = point_alpha, size = point_size,
                                                   ggplot2::aes(
                                                                text = paste0(
                                                                    "Time: ", time_k,
                                                                    "<br>Node: ", node,
                                                                    "<br>Outcome: ", outcome,
                                                                    "<br>Term: ", terms,
                                                                    "<br>Beta: ", signif(beta, 4)
                                                                ))))
    gg <- gg+ ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::labs(
                     x = "Models ordered by time and node group",
                     y = "Beta coefficient",
                     title = paste0("Regression coefficients of learner ",x$learner$name," across nuisance parameter models"),
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
    if (requireNamespace("plotly", quietly = TRUE)) {
        print(plotly::ggplotly(gg,tooltip = "text"))
    }else{
        print(gg)
    }
    invisible(gg)
}

######################################################################
### plot_nuisance.R ends here
