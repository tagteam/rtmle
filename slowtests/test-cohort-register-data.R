### test-cohort-register-data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: maj  1 2026 (09:51) 
## Version: 
## Last-Updated: maj  3 2026 (07:25) 
##           By: Thomas Alexander Gerds
##     Update #: 29
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(rtmle)
cohort <- simulate_cohort(
    n = 307,
    seed = 347,
    max_follow = 20,
    baseline_variables = list(age = "normal",sex = "binomial","SBP" = "normal"),
    intermediate_events = list(stroke = "Weibull",bleeding = "Weibull"),
    visit_events = list("A" = "binomial","B" = "binomial"),
    visit_measurements = list("changeSBP" = "normal"),
    visit_schedule = list(mean = 1,sd = 0.1,skip = 0,minimum_time_between_visits = 1),
    baseline_visit = list(A = "binomial",B = "binomial"),
    absorbing_events = list(death = "Weibull",dropout = "Weibull"),
    parameter_values = list(
        intercept_sex=qlogis(0.4),
        intercept_SBP=150,
        intercept_age=50,
        intercept_A=qlogis(0.3),
        intercept_B=qlogis(0.15),
        var_age=15,
        var_SBP=20,
        scale_death=0.01,
        scale_stroke=0.01,
        scale_bleeding=0.01,
        scale_dropout=0.01,
        effect_auto_A_A = 5,
        effect_auto_B_B = 4,
        effect_A_death = 0,
        effect_B_death = -0.8,
        effect_A_bleeding = 1.5,
        effect_B_bleeding = 0.5,
        effect_A_stroke = -0.5,
        effect_B_stroke = 0,
        effect_stroke_A = -2,
        effect_stroke_B = -1,
        effect_stroke_death = 2
    )
)
ld <- register_format(cohort)
x <- rtmle_init(time_grid = seq(0,20,4),
                name_id = "id",
                name_outcome = "stroke",
                name_competing = "death",
                name_censoring = "dropout",
                censored_label = "censored")
x <- add_long_data(x,
                   outcome_data=ld$timevar_data$stroke[!duplicated(id)],
                   censored_data=ld$timevar_data$dropout,
                   competing_data=ld$timevar_data$death,
                   timevar_data=ld$timevar_data[c("bleeding","changeSBP","A","B")])
x <- add_baseline_data(x,data=ld$baseline_data)
x <- long_to_wide(x,start_followup_date=0)
x <- prepare_rtmle_data(x)
x <- protocol(x,name = "Always_A",
              intervention = data.frame(time=x$intervention_nodes,
                                        "A" = factor("1",levels = c("0","1"))))
x <- protocol(x,name = "Never_A",
              intervention = data.frame(time=x$intervention_nodes,
                                        "A" = factor("0",levels = c("0","1"))))
x <- protocol(x,name = "Always_B",
              intervention = data.frame(time=x$intervention_nodes,
                                        "B" = factor("1",levels = c("0","1"))))
x <- target(x,name = "Treatment effect",
            protocols = c("Always_A","Always_B","Never_A"))
plot_adherence(x)
x <- model_formula(x)
x <- run_rtmle(x,time_horizon = 1:4)
x
######################################################################
### test-cohort-register-data.R ends here
