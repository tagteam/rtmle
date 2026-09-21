### test-dynamic-regimes.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: sep 12 2026 (07:08) 
## Version: 
## Last-Updated: sep 12 2026 (08:39) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# Load the checkout when this script is run from the source tree.  A plain
# `library(rtmle)` may otherwise attach an older installed copy with a stale
# API.  Fall back to the installed package when the script is copied outside
# a checkout.
source_root <- if (file.exists(file.path(getwd(), "DESCRIPTION"))) {
    normalizePath(getwd())
} else if (file.exists(file.path(getwd(), "..", "DESCRIPTION"))) {
    normalizePath(file.path(getwd(), ".."))
} else {
    NULL
}
if (!is.null(source_root) && requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(source_root, quiet = TRUE)
} else {
    library(rtmle)
}
data(simulated_cohort, package = "rtmle")
ld <- register_format(simulated_cohort)
y <- rtmle_init(time_grid = seq(0, 20, 4), name_id = "id",
                name_outcome = "stroke", name_competing = "death",
                name_censoring = "dropout", censored_label = "censored")
y <- add_long_data(
    y,
    outcome_data = ld$timevar_data$stroke[!duplicated(id)],
    censored_data = ld$timevar_data$dropout,
    competing_data = ld$timevar_data$death,
    timevar_data = ld$timevar_data[c("bleeding", "changeSBP", "A", "B")]
)
y <- add_baseline_data(y, data = ld$baseline_data)
y <- discretize_data(y, start_followup_date = 0)

# A dynamic intervention function returns the data after the intervention has
# been applied. The treatment propensity model describes the observed
# probability of following that history-dependent rule. Treatment at node k
# is based only on information available before that decision; the rule below
# uses the bleeding history through node k - 1.

has_bled_by <- function(data, node) {
    history <- intersect(paste0("bleeding_", 0:node), names(data))
    if (length(history) == 0L) {
        return(rep(FALSE, NROW(data)))
    }
    rowSums(do.call(cbind, lapply(history, function(v) {
        data[[v]] %in% c(1, "1")
    }))) > 0
}

stop_after_bleeding <- function(data, intervention_table, time_node) {
    intervened_data <- intervene(data, intervention_table, time_node)
    for (time in unique(intervention_table$time_node)) {
        treatment_node <- paste0("A_", time)
        has_bled <- has_bled_by(data, max(0, time - 1))
        if (is.factor(intervened_data[[treatment_node]])) {
            intervened_data[[treatment_node]][has_bled] <- factor(
                "0", levels = levels(intervened_data[[treatment_node]])
            )
        } else if (is.integer(intervened_data[[treatment_node]])) {
            intervened_data[[treatment_node]][has_bled] <- 0L
        } else {
            intervened_data[[treatment_node]][has_bled] <- 0
        }
    }
    intervened_data
}

stop_after_bleeding_A_and_B <- function(data, intervention_table, time_node) {
    intervened_data <- intervene(data, intervention_table, time_node)
    for (time in unique(intervention_table$time_node)) {
        has_bled <- has_bled_by(data, max(0, time - 1))
        for (treatment in c("A", "B")) {
            treatment_node <- paste0(treatment, "_", time)
            if (treatment_node %in% names(intervened_data)) {
                if (is.factor(intervened_data[[treatment_node]])) {
                    intervened_data[[treatment_node]][has_bled] <- factor(
                        "0", levels = levels(intervened_data[[treatment_node]])
                    )
                } else if (is.integer(intervened_data[[treatment_node]])) {
                    intervened_data[[treatment_node]][has_bled] <- 0L
                } else {
                    intervened_data[[treatment_node]][has_bled] <- 0
                }
            }
        }
    }
    intervened_data
}

# (a) Pooled adherence model: the observed treatment/intervention match is
# modelled using the ordinary pre-decision history.
y <- regime(
    y,
    name = "A_until_bleeding",
    intervention = data.frame(
        time_node = y$intervention_nodes,
        A = factor("1", levels = c("0", "1"))
    ),
    intervene_function = stop_after_bleeding,
    adherence_model_strata = NULL,
    verbose = FALSE
)

# (b) The same dynamic rule can be applied to two treatments. Factorization
# is a separate regime-level choice, not part of the adherence stratification.
y <- regime(
    y,
    name = "A1_B0_until_bleeding",
    intervention = data.frame(
        time_node = y$intervention_nodes,
        A = factor("1", levels = c("0", "1")),
        B = factor("0", levels = c("0", "1"))
    ),
    intervene_function = stop_after_bleeding_A_and_B,
    multiple_treatment_factorization = "sequential",
    adherence_model_strata = NULL,
    verbose = FALSE
)

# (c) Separate adherence models for patients with and without the latest
# pre-decision bleeding history. The unsuffixed stratum name is resolved to
# bleeding_(k - 1); precompute a cumulative history column if the regime
# concerns ever having bled rather than bleeding in the previous interval.
y <- regime(
    y,
    name = "A_until_bleeding_stratified",
    intervention = data.frame(
        time_node = y$intervention_nodes,
        A = factor("1", levels = c("0", "1"))
    ),
    intervene_function = stop_after_bleeding,
    adherence_model_strata = "bleeding",
    verbose = FALSE
)

# This regime deliberately keeps A = 1 after bleeding.  It is therefore a
# static Always-A regime and needs neither a custom intervention function nor
# bleeding in its propensity formula; the default intervene() function already
# applies A = 1 at every node.
y <- regime(
    y,
    name = "Always_A_even_after_bleeding",
    intervention = data.frame(
        time_node = y$intervention_nodes,
        A = factor("1", levels = c("0", "1"))
    ),
    verbose = FALSE
)

y <- prepare_rtmle_data(y)
y <- target(y, name = "Stroke_risk", estimator = "tmle",
            regimes = c("A_until_bleeding",
                          "A1_B0_until_bleeding",
                          "A_until_bleeding_stratified",
                          "Always_A_even_after_bleeding"))
y <- model_formula(y, verbose = FALSE)
y <- run_rtmle(y, learner = "learn_glm", time_horizon = 2,
               refit = TRUE, verbose = FALSE)
y


######################################################################
### test-dynamic-regimes.R ends here
