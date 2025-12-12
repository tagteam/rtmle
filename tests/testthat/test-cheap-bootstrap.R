### test-cheap-bootstrap.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr  9 2025 (10:02) 
## Version: 
## Last-Updated: Jun 17 2025 (07:28) 
##           By: Thomas Alexander Gerds
##     Update #: 17
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
library(data.table)
library(rtmle)
test_that("Cheap bootstrap confidence intervals",{
    tau <- 2
    set.seed(37)
    ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x <- add_long_data(x,
                    outcome_data=ld$outcome_data,
                    censored_data=ld$censored_data,
                    competing_data=ld$competing_data,
                    timevar_data=ld$timevar_data)
    x <- add_baseline_data(x,data=ld$baseline_data)
    x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
    x <- protocol(x,name = "Always_A",intervention = data.frame("A" = factor("1",levels = c("0","1"))),verbose = FALSE)
    x <- protocol(x,name = "Never_A",intervention = data.frame("A" = factor("0",levels = c("0","1"))),verbose = FALSE)
    x <- prepare_data(x)
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = c("Always_A","Never_A"))
    x <- model_formula(x)
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1,verbose = FALSE)
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau,verbose = FALSE)
    x <- cheap_bootstrap(x,B = 2,M = 71)
    a = x$estimate$Main_analysis[,.(Time_horizon,Protocol,Bootstrap_lower,Bootstrap_upper)]
    b = x$estimate$Cheap_bootstrap$Main_analysis[B == 2][,.(Time_horizon,Protocol,Bootstrap_lower,Bootstrap_upper)]
    setkey(a,Time_horizon,Protocol)
    setkey(b,Time_horizon,Protocol)
    expect_equal(a,b)
})


######################################################################
### test-cheap-bootstrap.R ends here
