### test-treat-unless.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: sep 20 2026 (07:33) 
## Version: 
## Last-Updated: sep 20 2026 (09:32) 
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
library(rtmle)
library(data.table)
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
    effect_GLP1_hypoglycemia = .3,
    effect_SGLT2_hypoglycemia = .2,
    effect_hypoglycemia_SGLT2 = -10,
    effect_hypoglycemia_GLP1 = -10,
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
    visit_events = list(
        GLP1 = "binomial", SGLT2 = "binomial"
    ),
    visit_measurements = list(deltaHbA1c = "normal"),
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
rf <- register_format(diabetes_cohort,treatment_variables = c("GLP1","SGLT2"))
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
    timevar_data = rf$timevar_data[c("hypoglycemia","af","deltaHbA1c",
                                     "GLP1","SGLT2")]
)
x <- discretize_data(
    x,
    # for long_data date formats provide a person-specific start date 
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

# The hand-written version below is kept as a reference implementation for
# comparison with treat_unless(). Define it before regime() because regime
# evaluates the intervention function while the regime is being registered.
stop_after_hypo <- function(data, intervention_table, time_node){
    # apply the underlying static regime
    # which sets treatment to 1
    out <- intervene(data, intervention_table, time_node)
    for (i in which(!is.na(intervention_table$value) &
                    intervention_table$time_node <= time_node &
                    grepl("^GLP1_", intervention_table$variable))) {
        v <- intervention_table$variable[i]
        # find history of the hypoglycemia process, including the current node
        h <- intersect(paste0("hypoglycemia_",
                              0:intervention_table$time_node[i]),
                       names(data))
        # change the value to zero for people
        # with a history of hypoglycemia
        if (length(h))
            out[[v]][rowSums(do.call(cbind, lapply(h, function(z)
                data[[z]] %in% c(1, "1")))) > 0] <- 0
    }
    out
}
x <- regime(
    x, name = "treat_until_hypo",
    intervention = data.frame(
        time_node=x$intervention_nodes,
        "GLP1" = factor(rep("1",10),levels = c("0","1"))
    ), 
    intervene_function = treat_unless(
        contra_indication = "hypoglycemia",
        action = 0,
        lookback_window = Inf
    )
)
x <- regime(
    x, name = "treat-until-hypo",
    intervention = data.frame(
        time_node=x$intervention_nodes,
        "GLP1" = factor(rep("1",10),levels = c("0","1"))
    ), 
    intervene_function = "stop_after_hypo"
)

stopifnot(identical(
    x$regimes[["treat_until_hypo"]]$intervention_match,
    x$regimes[["treat-until-hypo"]]$intervention_match
))
a <- adherence(x)
summary_adherence(x)
## plot_adherence(x)

######################################################################
### test-treat-unless.R ends here
