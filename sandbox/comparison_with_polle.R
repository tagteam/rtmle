### comparison_with_polle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: sep 21 2026 (15:42) 
## Version: 
## Last-Updated: sep 22 2026 (09:56) 
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
library(rtmle)
library(data.table)
library(polle)
diabetes_parameters <- list(
  intercept_age = 60, var_age = 10, intercept_sex = 1,
  intercept_HbA1c = 58, var_HbA1c = 8,
  intercept_GLP1 = 1, intercept_SGLT2 = 1,
  intercept_deltaHbA1c = 0, var_deltaHbA1c = 2,
  scale_MACE = .0015, scale_death = .0007, scale_dropout = .0012,
  scale_af = 0.01, scale_hypoglycemia = .0015, 
  effect_age_MACE = .01, effect_HbA1c_MACE = .02,
  effect_GLP1_MACE = -.3, effect_SGLT2_MACE = -.3,
  effect_GLP1_death = -.1, effect_SGLT2_death = -.1,
  effect_GLP1_dropout = .1, effect_SGLT2_dropout = .1,
  effect_GLP1_hypoglycemia = .3, effect_hypoglycemia_GLP1 = -10,
  effect_SGLT2_hypoglycemia = .2, effect_hypoglycemia_SGLT2 = -3,
  effect_GLP1_deltaHbA1c = -1.5,
  effect_SGLT2_deltaHbA1c = -1.2,
  effect_age_GLP1 = .005, effect_HbA1c_GLP1 = .01,
  effect_age_SGLT2 = .005, effect_HbA1c_SGLT2 = .01
)
diabetes_cohort <- simulate_cohort(
  n = 1000, seed = 20260915,
  # clock is in months
  max_follow = 60,
  # baseline
  baseline_variables = list(age = "normal", sex = "binomial", HbA1c = "normal"),
  baseline_visit = list(
    GLP1 = "binomial", SGLT2 = "binomial"
  ),
  # hooks make the function flexible
  post_baseline_visit_hook = function(X) {
    X[, `:=`(GLP1 = rbinom(.N, 1, .45),SGLT2 = rbinom(.N, 1, .30))]
  },
  # visit schedule
  visit_schedule = list(
    mean = 6, sd = .2, skip = 0,
    minimum_time_between_visits = 1
  ),
  # blood test
  visit_measurements = list(deltaHbA1c = "normal"),
  # treatment
  visit_events = list(
    GLP1 = "binomial", SGLT2 = "binomial"
  ),
  # comorbidity and side effects
  intermediate_events = list(
    af = "Weibull",
    hypoglycemia = "Weibull"
  ),
  # end-of-story
  absorbing_events = list(
    MACE = "Weibull",
    death = "Weibull",
    dropout = "Weibull"
  ),
  # regression parameters
  parameter_values = diabetes_parameters
)
rf <- register_format(
  diabetes_cohort,
  treatment_variables = c("GLP1","SGLT2")
)
x <- rtmle_init(
  time_unit = "Months since treatment start",
  time_grid = seq(0, 60, 6),
  name_id = "id",
  name_outcome = "MACE",
  name_competing = "death",
  name_censoring = "dropout",
  censored_levels = c("uncensored","censored"),
  censored_label = "censored",
  # hyperparameters
  minority_threshold = 8, 
  weight_truncation = c(0, 1),
  prediction_range = c(0.0001, 0.9999)
)
# add data to the object
x <- add_baseline_data(
  x,
  data = rf$baseline_data
)
x <- add_long_data(
  x,
  outcome_data = rf$timevar_data$MACE,
  censored_data = rf$timevar_data$dropout,
  competing_data = rf$timevar_data$death,
  timevar_data = rf$timevar_data[c("hypoglycemia","af","deltaHbA1c","GLP1","SGLT2")]
)
x <- discretize_data(
  x,
  # date format: provide a person-specific start date 
  start_followup_date = 0,
  hypoglycemia = list(method = "event"),
  af = list(method = "event"),
  deltaHbA1c = list(method = "locf"),
  GLP1 = list(
    variable = "GLP1",
    method = "exposure_percent",
    threshold = 0.8
  ),
  SGLT2 = list(
    variable = "SGLT2",
    method = "exposure_percent",
    threshold = 0.8
  )
)
x <- prepare_rtmle_data(x)
x <- protocol(
  x, name = "always_GLP1_never_SGLT2",
  intervention = data.frame(
    time_node=x$intervention_nodes,
    "GLP1" = factor(rep("1",10),levels = c("0","1")),
    "SGLT2" = factor(rep("0",10),levels = c("0","1"))
  )
)
# short form
x <- protocol(
  x, name = "always_SGLT2_never_GLP1", 
  intervention = c(1,0),
  treatment_variables = c("SGLT2", "GLP1")
)
# Partial control: NA means 'no intervention'
x <- protocol(
  x, name = "limit_GLP1_allow_SGLT2",
  intervention = data.frame(
    time_node=x$intervention_nodes,
    "GLP1" = factor(c(rep("1",5),rep(NA,5)),levels = c("0","1")),
    "SGLT2" = factor(c(0,0,rep(NA,8)),levels = c("0","1"))
  ),
  # static rule
  intervene_function = "intervene"
)
x <- protocol(
  x, name = "treat_until_hypo",
  intervention = data.frame(
    time_node=x$intervention_nodes,
    "GLP1" = factor(rep("1",10),levels = c("0","1")),
    "SGLT2" = factor(rep("0",10),levels = c("0","1"))
  ),
  intervene_function = treat_unless(
    contra_indication = "hypoglycemia",
    action = c("0",NA)
  )
)
x <- target(x, name = "ateMACE", estimator = "tmle",
            regimes = c("always_GLP1_never_SGLT2",
                          "always_SGLT2_never_GLP1"))
x <- target(x, name = "suppl_MACE_risk", estimator = "tmle",
            regimes = c("limit_GLP1_allow_SGLT2",
                          "treat_until_hypo"))
x <- model_formula(x,exclusion_rules = list("GLP1_*" = "SGLT2_0"),Markov = c("af","hypoglycemia","deltaHbA1c"))
x <- run_rtmle(x,refit = FALSE,time_horizon = 1:8,learner = "learn_glmnet")
x


######################################################################
### The same target parameters with polle
######################################################################
#
# polle works with a finite-stage policy value.  The conversion below uses
# two action rows for every rtmle intervention node and a terminal utility
#
#     U = 1{first terminal event is MACE}.
#
# Death and administrative end of follow-up are terminal events with U = 0;
# dropout is event = 2 (right censoring).  Thus E[U] is the cumulative MACE
# incidence in the presence of death as a competing terminal event.
#
# There are two details worth making explicit.
#
# 1. polle has one categorical action per stage.  The two rtmle treatments are
#    therefore represented by two consecutive action stages at each node:
#    GLP1 first, followed by SGLT2.  This is an ordered factorization of the
#    same joint treatment assignment, not a change in the intervention.
# 2. A treatment component left natural is represented by the one-level action
#    "natural" and its observed value is carried in the state.  This preserves
#    the meaning of NA in the rtmle intervention.  It also avoids multiplying
#    four-level joint-treatment propensities when two binary treatments are
#    both controlled.
#
# Two dummy natural actions after the requested horizon keep all horizon-end
# censoring rows before polle's final utility stage.  This avoids asking polle
# to impute the final utility and lets its right-censoring score use the
# observed terminal utility whenever it is available.

make_polle_stage_data <- function(prepared_data,
                                  horizon,
                                  regime_name,
                                  baseline_names = c("age", "sex", "HbA1c")) {
    stopifnot(length(horizon) == 1L, horizon >= 1L,
              horizon == as.integer(horizon))
    horizon <- as.integer(horizon)
    d <- data.table::copy(data.table::as.data.table(prepared_data))

    value01 <- function(name, i, node) {
        as.numeric(as.character(d[[paste0(name, "_", node)]][i]))
    }

    first_event <- function(i) {
        for (j in seq_len(horizon)) {
            if (isTRUE(d[[paste0("MACE_", j)]][i] == 1)) {
                return(list(type = "MACE", time = j))
            }
            if (isTRUE(d[[paste0("death_", j)]][i] == 1)) {
                return(list(type = "death", time = j))
            }
            if (isTRUE(as.character(
                d[[paste0("dropout_", j)]][i]
            ) == "censored")) {
                return(list(type = "dropout", time = j))
            }
        }
        list(type = NA_character_, time = NA_integer_)
    }

    has_hypoglycemia <- function(i, node) {
        if (node < 0L) return(FALSE)
        any(vapply(seq_len(node + 1L), function(j) {
            value01("hypoglycemia", i, j - 1L) != 0
        }, logical(1)))
    }

    controlled_components <- function(i, node) {
        if (regime_name %in% c("always_GLP1_never_SGLT2",
                               "always_SGLT2_never_GLP1")) {
            return(c(GLP1 = TRUE, SGLT2 = TRUE))
        }
        if (regime_name == "limit_GLP1_allow_SGLT2") {
            return(c(GLP1 = node < 5L, SGLT2 = node < 2L))
        }
        if (regime_name == "treat_until_hypo") {
            return(c(GLP1 = TRUE,
                     SGLT2 = !has_hypoglycemia(i, node)))
        }
        stop("Unknown regime: ", regime_name)
    }

    state_at_node <- function(i, node, component, controlled,
                              previous_GLP1, previous_SGLT2) {
        glp1 <- value01("GLP1", i, node)
        sglt2 <- value01("SGLT2", i, node)
        data.table::data.table(
            hypoglycemia = value01("hypoglycemia", i, node),
            af = value01("af", i, node),
            deltaHbA1c = value01("deltaHbA1c", i, node),
            hypo_ever = as.integer(has_hypoglycemia(i, node)),
            natural_GLP1 = if (controlled[["GLP1"]]) 0 else glp1,
            natural_SGLT2 = if (controlled[["SGLT2"]]) 0 else sglt2,
            lag_GLP1 = if (component == "GLP1") {
                previous_GLP1
            } else {
                glp1
            },
            lag_SGLT2 = previous_SGLT2
        )
    }

    one_subject <- function(i) {
        terminal <- first_event(i)
        last_node_count <- if (is.na(terminal$time)) horizon else terminal$time
        rows <- vector("list", 0L)

        for (node in 0:(last_node_count - 1L)) {
            controlled <- controlled_components(i, node)
            glp1 <- value01("GLP1", i, node)
            sglt2 <- value01("SGLT2", i, node)
            previous_GLP1 <- if (node == 0L) 0 else value01(
                "GLP1", i, node - 1L
            )
            previous_SGLT2 <- if (node == 0L) 0 else value01(
                "SGLT2", i, node - 1L
            )
            rows[[length(rows) + 1L]] <- c(
                data.table::data.table(
                    id = d$id[i], stage = 2L * node + 1L,
                    event = 0L,
                    A = if (controlled[["GLP1"]]) as.character(glp1)
                        else "natural",
                    U = 0, time = node, time2 = node + 1L
                ),
                state_at_node(
                    i, node, "GLP1", controlled,
                    previous_GLP1, previous_SGLT2
                )
            )
            rows[[length(rows) + 1L]] <- c(
                data.table::data.table(
                    id = d$id[i], stage = 2L * node + 2L,
                    event = 0L,
                    A = if (controlled[["SGLT2"]]) as.character(sglt2)
                        else "natural",
                    U = 0, time = node, time2 = node + 1L
                ),
                state_at_node(
                    i, node, "SGLT2", controlled,
                    previous_GLP1, previous_SGLT2
                )
            )
        }

        if (is.na(terminal$time)) {
            # Two natural bookkeeping actions keep horizon-end censoring rows
            # before polle's final utility stage.
            node <- horizon - 1L
            previous_GLP1 <- if (node == 0L) 0 else value01(
                "GLP1", i, node - 1L
            )
            previous_SGLT2 <- if (node == 0L) 0 else value01(
                "SGLT2", i, node - 1L
            )
            natural <- c(GLP1 = FALSE, SGLT2 = FALSE)
            rows[[length(rows) + 1L]] <- c(
                data.table::data.table(
                    id = d$id[i], stage = 2L * horizon + 1L,
                    event = 0L, A = "natural", U = 0,
                    time = horizon, time2 = horizon + 1L
                ),
                state_at_node(
                    i, node, "GLP1", natural,
                    previous_GLP1, previous_SGLT2
                )
            )
            rows[[length(rows) + 1L]] <- c(
                data.table::data.table(
                    id = d$id[i], stage = 2L * horizon + 2L,
                    event = 0L, A = "natural", U = 0,
                    time = horizon, time2 = horizon + 1L
                ),
                state_at_node(
                    i, node, "SGLT2", natural,
                    previous_GLP1, previous_SGLT2
                )
            )
            rows[[length(rows) + 1L]] <- c(
                data.table::data.table(
                    id = d$id[i], stage = 2L * horizon + 3L,
                    event = 1L, A = NA_character_, U = 0,
                    time = horizon + 1L, time2 = horizon + 2L
                ),
                state_at_node(
                    i, node, "SGLT2", natural,
                    previous_GLP1, previous_SGLT2
                )
            )
        } else {
            node <- last_node_count - 1L
            controlled <- controlled_components(i, node)
            previous_GLP1 <- if (node == 0L) 0 else value01(
                "GLP1", i, node - 1L
            )
            previous_SGLT2 <- if (node == 0L) 0 else value01(
                "SGLT2", i, node - 1L
            )
            endpoint_event <- if (terminal$type == "dropout") 2L else 1L
            endpoint_utility <- if (terminal$type == "MACE") 1 else 0
            rows[[length(rows) + 1L]] <- c(
                data.table::data.table(
                    id = d$id[i], stage = 2L * last_node_count + 1L,
                    event = endpoint_event, A = NA_character_,
                    # For censoring before K+1 this placeholder is ignored by
                    # get_utility(), which marks the subject as censored.
                    U = endpoint_utility,
                    time = last_node_count,
                    time2 = last_node_count + 1L
                ),
                state_at_node(
                    i, node, "SGLT2", controlled,
                    previous_GLP1, previous_SGLT2
                )
            )
        }
        data.table::rbindlist(rows, fill = TRUE)
    }

    list(
        stage_data = data.table::rbindlist(
            lapply(seq_len(nrow(d)), one_subject),
            fill = TRUE
        ),
        baseline_data = d[, c("id", baseline_names), with = FALSE],
        action_set = c("0", "1", "natural")
    )
}

make_polle_policy <- function(regime_name, horizon) {
    horizon <- as.integer(horizon)
    K <- 2L * horizon + 2L
    if (regime_name == "always_GLP1_never_SGLT2") {
        actions <- c(rep(c("1", "0"), horizon), "natural", "natural")
        return(polle::policy_def(actions, name = regime_name))
    }
    if (regime_name == "always_SGLT2_never_GLP1") {
        actions <- c(rep(c("0", "1"), horizon), "natural", "natural")
        return(polle::policy_def(actions, name = regime_name))
    }
    if (regime_name == "limit_GLP1_allow_SGLT2") {
        actions <- unlist(lapply(0:(horizon - 1L), function(node) {
            if (node < 2L) c("1", "0")
            else if (node < 5L) c("1", "natural")
            else c("natural", "natural")
        }))
        actions <- c(actions, "natural", "natural")
        return(polle::policy_def(actions, name = regime_name))
    }
    if (regime_name == "treat_until_hypo") {
        rule <- function(hypo_ever, stage, ...) {
            node <- (stage - 1L) %/% 2L
            ifelse(stage > 2L * horizon, "natural",
                   ifelse(stage %% 2L == 1L,
                          ifelse(hypo_ever == 0, "1", "0"),
                          ifelse(hypo_ever == 0, "0", "natural")))
        }
        return(polle::policy_def(
            rep(list(rule), K),
            stage_number = TRUE,
            name = regime_name
        ))
    }
    stop("Unknown regime: ", regime_name)
}

make_polle_models <- function(policy_data) {
    action_counts <- vapply(
        policy_data$stage_action_sets,
        length,
        integer(1)
    )
    g_models <- lapply(action_counts, function(n_actions) {
        if (n_actions > 1L) {
            polle::g_glm(
                formula = ~hypo_ever + age + sex + HbA1c +
                    lag_GLP1 + lag_SGLT2
                ## num.trees = 200,
                ## min.node.size = 30
            )
        } else {
            polle::g_empir(formula = ~1)
        }
    })
    q_models <- lapply(action_counts, function(n_actions) {
        if (n_actions > 1L) {
            polle::q_glm(formula = ~A + .)
        } else {
            polle::q_glm(formula = ~1)
        }
    })
    c_models <- replicate(
        length(action_counts) + 1L,
        polle::g_empir(formula = ~1),
        simplify = FALSE
    )
    list(g = g_models, q = q_models, c = c_models)
}

polle_value <- function(policy_data, policy, models, name) {
    polle::policy_eval(
        policy_data = policy_data,
        policy = policy,
        g_models = models$g,
        q_models = models$q,
        c_models = models$c,
        type = "dr",
        name = name
    )
}

polle_value_summary <- function(pe, regime_name, horizon) {
    ic <- as.numeric(polle::IC(pe))
    estimate <- as.numeric(stats::coef(pe))
    data.table::data.table(
        Regime = regime_name,
        Time_horizon = horizon,
        Estimate = estimate,
        Standard_error = sqrt(mean((ic - mean(ic))^2) / length(ic))
    )
}

polle_regimes <- c(
    "always_GLP1_never_SGLT2",
    "always_SGLT2_never_GLP1",
    "limit_GLP1_allow_SGLT2",
    "treat_until_hypo"
)

# Make the random-forest propensity fits reproducible.
set.seed(20260915)
polle_fits <- list()
polle_estimates <- data.table::rbindlist(lapply(seq_len(8), function(horizon) {
    rows <- lapply(polle_regimes, function(regime_name) {
        formatted <- make_polle_stage_data(
            x$prepared_data,
            horizon = horizon,
            regime_name = regime_name
        )
        pd <- polle::policy_data(
            data = formatted$stage_data,
            baseline_data = formatted$baseline_data,
            type = "long",
            id = "id",
            stage = "stage",
            event = "event",
            action = "A",
            utility = "U",
            time = "time",
            time2 = "time2",
            action_set = formatted$action_set
        )
        stopifnot(polle::get_K(pd) == 2L * horizon + 2L)
        models <- make_polle_models(pd)
        pe <- polle_value(
            pd,
            make_polle_policy(regime_name, horizon),
            models,
            regime_name
        )
        polle_fits[[paste(regime_name, horizon, sep = "__")]] <<- list(
            policy_data = pd,
            policy = make_polle_policy(regime_name, horizon),
            fit = pe
        )
        polle_value_summary(pe, regime_name, horizon)
    })
    data.table::rbindlist(rows)
}), fill = TRUE)

polle_contrasts <- data.table::rbindlist(lapply(seq_len(8), function(horizon) {
    value_for <- function(regime_name) {
        polle_fits[[paste(regime_name, horizon, sep = "__")]]$fit
    }
    contrast_for <- function(regime, reference, target_parameter) {
        p <- value_for(regime)
        r <- value_for(reference)
        zp <- as.numeric(polle::IC(p)) + as.numeric(stats::coef(p))
        zr <- as.numeric(polle::IC(r)) + as.numeric(stats::coef(r))
        mp <- mean(zp)
        mr <- mean(zr)
        if (target_parameter == "Risk_difference") {
            z <- zp - zr
            estimate <- mp - mr
        } else {
            z <- (zp - mp) / mr - mp * (zr - mr) / (mr^2)
            estimate <- mp / mr
        }
        data.table::data.table(
            Regime = regime,
            Reference = reference,
            Target_parameter = target_parameter,
            Time_horizon = horizon,
            Estimate = estimate,
            Standard_error = sqrt(mean((z - mean(z))^2) / length(z))
        )
    }
    data.table::rbindlist(list(
        contrast_for("always_SGLT2_never_GLP1",
                     "always_GLP1_never_SGLT2", "Risk_difference"),
        contrast_for("always_SGLT2_never_GLP1",
                     "always_GLP1_never_SGLT2", "Risk_ratio"),
        contrast_for("treat_until_hypo",
                     "limit_GLP1_allow_SGLT2", "Risk_difference"),
        contrast_for("treat_until_hypo",
                     "limit_GLP1_allow_SGLT2", "Risk_ratio")
    ))
}), fill = TRUE)

# Put the polle output beside the corresponding rtmle output.  rtmle stores
# contrasts without a separate Reference column, so add the reference implied
# by the two target definitions before joining the results.
rtmle_comparison <- data.table::copy(x$estimate$Main_analysis)
rtmle_comparison[, `:=`(
    Target = as.character(Target),
    Regime = as.character(Regime),
    Target_parameter = as.character(Target_parameter)
)]
rtmle_comparison[, Reference := NA_character_]
rtmle_comparison[
    Target_parameter %in% c("Risk_difference", "Risk_ratio") &
        Target == "ateMACE",
    Reference := "always_GLP1_never_SGLT2"
]
rtmle_comparison[
    Target_parameter %in% c("Risk_difference", "Risk_ratio") &
        Target == "suppl_MACE_risk",
    Reference := "limit_GLP1_allow_SGLT2"
]
rtmle_comparison <- rtmle_comparison[, .(
    Target, Regime, Reference, Target_parameter, Time_horizon,
    rtmle_Estimate = Estimate,
    rtmle_Standard_error = Standard_error
)]

polle_comparison <- data.table::rbindlist(list(
    polle_estimates[, .(
        Target = ifelse(
            as.character(Regime) %in% c(
                "always_GLP1_never_SGLT2",
                "always_SGLT2_never_GLP1"
            ),
            "ateMACE", "suppl_MACE_risk"
        ),
        Regime = as.character(Regime),
        Reference = NA_character_,
        Target_parameter = "Risk",
        Time_horizon,
        polle_Estimate = Estimate,
        polle_Standard_error = Standard_error
    )],
    polle_contrasts[, .(
        Target = ifelse(
            as.character(Regime) == "always_SGLT2_never_GLP1",
            "ateMACE", "suppl_MACE_risk"
        ),
        Regime = as.character(Regime),
        Reference = as.character(Reference),
        Target_parameter = as.character(Target_parameter),
        Time_horizon,
        polle_Estimate = Estimate,
        polle_Standard_error = Standard_error
    )]
), fill = TRUE)
comparison <- merge(
    rtmle_comparison,
    polle_comparison,
    by = c("Target", "Regime", "Reference", "Target_parameter",
           "Time_horizon"),
    all = TRUE
)
comparison[, Difference := polle_Estimate - rtmle_Estimate]

# Keep the comparison attached to the rtmle object for interactive inspection.
x$polle <- list(
    estimates = polle_estimates,
    contrasts = polle_contrasts,
    comparison = comparison,
    fits = polle_fits
)

polle_estimates
polle_contrasts
comparison


######################################################################
### comparison_with_polle.R ends here
