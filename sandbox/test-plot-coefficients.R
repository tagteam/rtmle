### test-plot-coefficients.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 23 2026 (10:42) 
## Version: 
## Last-Updated: mar 18 2026 (08:17) 
##           By: Thomas Alexander Gerds
##     Update #: 5
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# -------------------------
# Example usage:
#
# Single panel "manhattan style", all terms, colored by node group:

z <- run_analysis(tau = 1:11,intervals = c(0,30.45,3*30.45,seq(6*30.45,9*6*30.45,6*30.45)),outcome = "primary.outcome",dt_baseline = dt_baseline,discretized_treatment_regimes = discretized_treatment_regimes,discretized_timevar_covariates = discretized_timevar_covariates,discretized_outcomes = discretized_outcomes,names_baseline_covariates = names_baseline_covariates,timevarying_covariate_groups = timevarying_covariate_groups,markov_variables = NULL,no_lira = FALSE,learner = list(name = "ridge",alpha = 1,selector = "undersmooth",fun = "learn_glmnet"))
z <- run_analysis(tau = 1:11,intervals = c(0,30.45,3*30.45,seq(6*30.45,9*6*30.45,6*30.45)),outcome = "primary.outcome",dt_baseline = dt_baseline,discretized_treatment_regimes = discretized_treatment_regimes,discretized_timevar_covariates = discretized_timevar_covariates,discretized_outcomes = discretized_outcomes,names_baseline_covariates = names_baseline_covariates,timevarying_covariate_groups = timevarying_covariate_groups,markov_variables = NULL,no_lira = FALSE,learner = list(name = "lasso",alpha = 0,selector = "undersmooth",fun = "learn_glmnet"))

res <- plot_model_coefficients(
  z,
  plot_style = "manhattan",
  times = 0:10,
  node_order = c("protocol","censoring","outcome"),
  protocol_nodes = c("Placebo","Lira"),
  manhattan_color_by = "node_group",
  show_x_labels = FALSE
)
ggplotly(res$plot, tooltip = "text")

# res$plot
#
# If it's too dense, filter terms:
# res2 <- plot_model_betas(x, plot_style="manhattan", term="sexMale|age60-80", manhattan_color_by="time")
# res2$plot
p <- ggplotly(res$plot, tooltip = "text")
htmlwidgets::saveWidget(p, "~/tmp/betas.html", selfcontained = TRUE)
browseURL("~/tmp/betas.html")

library(ggplot2)
library(plotly)
library(htmlwidgets)

d <- data.frame(x = 1:10, y = rnorm(10))
g <- ggplot(d, aes(x, y, text = paste0("pt ", x, ": ", round(y, 3)))) + geom_point()



p <- ggplotly(g, tooltip = "text")
saveWidget(p, "~/tmp/sanity.html", selfcontained = TRUE)
browseURL("~/tmp/sanity.html")

# How many traces?
length(res$data)

# Show whether each trace has text + hoverinfo
info <- lapply(res$data, function(tr) {
  list(
    type = tr$type,
    mode = tr$mode,
    has_text = !is.null(tr$text),
    text_head = if (!is.null(tr$text)) head(tr$text, 2) else NULL,
    hoverinfo = tr$hoverinfo
  )
})
str(info, 2)


######################################################################
### test-plot-coefficients.R ends here
