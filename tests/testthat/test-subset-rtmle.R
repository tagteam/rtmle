### test-subset-rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr  8 2025 (14:47) 
## Version: 
## Last-Updated: maj  7 2026 (16:29) 
##           By: Thomas Alexander Gerds
##     Update #: 47
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
test_that("run rtmle on a subset",{
    set.seed(17)
    ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(time_grid = seq(0,1500,30.45*12),name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x <- add_long_data(x,outcome_data=ld$outcome_data,censored_data=ld$censored_data,competing_data=ld$competing_data,timevar_data=ld$timevar_data)
    x <- add_baseline_data(x,data=ld$baseline_data)
    x <- discretize_data(x,start_followup_date = 0)
    x <- regime(x,name = "Always_A",intervention = data.frame("time"=x$intervention_nodes,"A" = factor("1",levels = c("0","1"))),verbose = 0L)
    x <- regime(x,name = "Never_A",intervention = data.frame("time"=x$intervention_nodes,"A" = factor("0",levels = c("0","1"))),verbose = 0L)
    x <- prepare_rtmle_data(x)
    x <- target(x,name = "Outcome_risk",estimator = "tmle",regimes = c("Never_A","Always_A"))
    x <- target(x,name = "Test Outcome_risk",estimator = "tmle",regimes = c("Never_A"))
    x <- model_formula(x)
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:2,verbose = 0L)
    # subset to first 79
    ld1 <- lapply(ld,function(d){subset(d,d$id <= 79)})
    names(ld1) <- names(ld)
    ld1$timevar_data <- lapply(ld$timevar_data,function(d){subset(d,d$id <= 79)})
    names(ld1$timevar_data) <- names(ld$timevar_data)
    # now run again on subset
    x1 <- rtmle_init(time_grid = seq(0,1500,30.45*12),name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x1 <- add_long_data(x1,outcome_data=ld1$outcome_data,censored_data=ld1$censored_data,competing_data=ld1$competing_data,timevar_data=ld1$timevar_data)
    x1 <- add_baseline_data(x1,data = ld1$baseline_data)
    x1 <- discretize_data(x1,start_followup_date = 0)
    x1 <- regime(x1,name = "Always_A",intervention = data.frame("time"=x$intervention_nodes,"A" = factor("1",levels = c("0","1"))),verbose = 0L)
    x1 <- regime(x1,name = "Never_A",intervention = data.frame("time"=x$intervention_nodes,"A" = factor("0",levels = c("0","1"))),verbose = 0L)
    x1 <- prepare_rtmle_data(x1)
    x1 <- target(x1,name = "Outcome_risk",estimator = "tmle",regimes = c("Never_A","Always_A"))
    x1 <- target(x1,name = "Test Outcome_risk",estimator = "tmle",regimes = c("Never_A"))
    x1 <- model_formula(x1)
    x1 <- run_rtmle(x1,learner = "learn_glmnet",time_horizon = 2,verbose = 0L)
    x79 <- run_rtmle(subset(x, ids = 1:79), learner = "learn_glmnet", time_horizon = 2)
    expect_equal(x1$estimate$Main_analysis, x79$estimate$Main_analysis)
    expect_equal(x1$IC, x79$IC)
})

test_that("stratified analyses",
{
    set.seed(17)
    ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(time_grid = seq(0,1500,30.45*12),name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x <- add_long_data(x,outcome_data=ld$outcome_data,censored_data=ld$censored_data,competing_data=ld$competing_data,timevar_data=ld$timevar_data)
    x <- add_baseline_data(x,data=ld$baseline_data)
    x <- discretize_data(x,start_followup_date = 0)
    x <- regime(x,name = "Always_A",intervention = data.frame("time"=x$intervention_nodes,"A" = factor("1",levels = c("0","1"))),verbose = 0L)
    x <- regime(x,name = "Never_A",intervention = data.frame("time"=x$intervention_nodes,"A" = factor("0",levels = c("0","1"))),verbose = 0L)
    x <- prepare_rtmle_data(x)
    x <- target(x,name = "Outcome_risk",estimator = "tmle",regimes = c("Always_A","Never_A"))
    x <- model_formula(x)
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:2,verbose = 0L)
    strata <- list(Female = subset(x, sex == 0), Male = subset(x, sex == 1))
    fits <- lapply(strata, run_rtmle, learner = "learn_glmnet", time_horizon = 2)
    combined_summary <- data.table::rbindlist(lapply(fits, summary), idcol = "Sex")
    expect_setequal(unique(combined_summary$Sex), c("Female", "Male"))
    expect_false("Sex" %in% names(summary(x)))
})

######################################################################
### test-subset-rtmle.R ends here
