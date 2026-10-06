### test-regime.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: sep 23 2026 (08:47) 
## Version: 
## Last-Updated: sep 23 2026 (17:01) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
library(rtmle)
library(data.table)
test_that("interplay between regime and model_formula works with arguments adherence_model_strata and multiple_treatment_factorization",{
    data(simulated_cohort)
    rf <- register_format(simulated_cohort,treatment_variables = c("A","B"))
    x <- rtmle_init(time_unit = "Months since treatment start",time_grid = seq(0,18,3),name_id = "id",name_outcome = "death",name_competing = NULL,name_censoring = "dropout",censored_levels = c("uncensored","censored"),censored_label = "censored",minority_threshold = 8,weight_truncation = c(0, 1),prediction_range = c(0.0001, 0.9999))
    # providing a non-existing variable as strata should not yet break the code
    x <- protocol(
        x,
        name = "break",
        intervention = data.frame(time_node=x$intervention_nodes,"A" = factor(rep("1",6),levels = c("0","1")),"B" = factor(rep("0",6),levels = c("0","1"))),
        intervene_function = treat_unless(contra_indication = "bleeding",action = c("0",NA)),
        multiple_treatment_factorization = "independent",
        adherence_model_strata = "something"
    )
    # but model_formula should throw an error
    x <- model_formula(x,exclusion_rules = list("A_*" = "B_0"),Markov = c("bleeding","stroke","changeSBP"))
    x <- add_baseline_data(x,data = rf$baseline_data)
    x <- add_long_data(x,outcome_data = rf$timevar_data$death,censored_data = rf$timevar_data$dropout,competing_data = NULL,timevar_data = rf$timevar_data[c("stroke","bleeding","changeSBP","A","B")])
    x <- discretize_data(x,start_followup_date = 0,stroke = list(method = "event"),bleeding = list(method = "event"),changeSBP = list(method = "locf"),A = list(variable = "A",method = "exposure_percent",threshold = 0.8),B = list(variable = "B",method = "exposure_percent",threshold = 0.8))
    x <- prepare_rtmle_data(x)
    x <- protocol(
        x,
        name = "treat_until_bleeding_joint",
        intervention = data.frame(time_node=x$intervention_nodes,"A" = factor(rep("1",6),levels = c("0","1")),"B" = factor(rep("0",6),levels = c("0","1"))),
        intervene_function = treat_unless(contra_indication = "bleeding",action = c("0",NA)),
        multiple_treatment_factorization = "joint",
        adherence_model_strata = NULL
    )
    x <- protocol(
        x,
        name = "treat_until_bleeding_seq",
        intervention = data.frame(time_node=x$intervention_nodes,"A" = factor(rep("1",6),levels = c("0","1")),"B" = factor(rep("0",6),levels = c("0","1"))),
        intervene_function = treat_unless(contra_indication = "bleeding",action = c("0",NA)),
        multiple_treatment_factorization = "sequential",
        adherence_model_strata = "bleeding"
    )
    x <- protocol(
        x,
        name = "treat_until_bleeding_ind",
        intervention = data.frame(time_node=x$intervention_nodes,"A" = factor(rep("1",6),levels = c("0","1")),"B" = factor(rep("0",6),levels = c("0","1"))),
        intervene_function = treat_unless(contra_indication = "bleeding",action = c("0",NA)),
        multiple_treatment_factorization = "independent",
        adherence_model_strata = "bleeding"
    )
    # providing a non-existing variable as strata should now break the code
    expect_error(x <- protocol(
                     x,
                     name = "break",
                     intervention = data.frame(time_node=x$intervention_nodes,"A" = factor(rep("1",6),levels = c("0","1")),"B" = factor(rep("0",6),levels = c("0","1"))),
                     intervene_function = treat_unless(contra_indication = "bleeding",action = c("0",NA)),
                     multiple_treatment_factorization = "independent",
                     adherence_model_strata = "something"
                 ))
    x <- model_formula(x,exclusion_rules = list("A_*" = "B_0"),Markov = c("bleeding","stroke","changeSBP"))
    x$models$time_2$treat_until_bleeding_joint
    x$models$time_2$treat_until_bleeding_seq
    x$models$time_2$treat_until_bleeding_ind
})


######################################################################
### test-regime.R ends here
