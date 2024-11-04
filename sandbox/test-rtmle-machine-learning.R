### test-rtmle-machine-learning.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 28 2024 (10:19) 
## Version: 
## Last-Updated: Nov  2 2024 (07:42) 
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
library(targets)
library(prodlim)
set.seed(17)
ld <- simulate_long_data(n = 1113,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
protocol(x) <- list(name = "Always_A",treatment_variables = "A",intervention = 1)
prepare_data(x) <- list()
target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "Always_A")
system.time(x <- run_rtmle(x,refit = TRUE,time_horizon = 2,folds = 10,learner = c("learn_glm","learn_glm")))
summary(x)
system.time(x <- run_rtmle(x,refit = TRUE,time_horizon = 1:3,learner = "learn_glm"))
summary(x)
system.time(x <- run_rtmle(x,refit = TRUE,time_horizon = 1:3,learner = "learn_glmnet"))
summary(x)
system.time(x <- run_rtmle(x,refit = TRUE,time_horizon = 2,learner = "learn_ranger"))
summary(x)
system.time(x <- run_rtmle(x,refit = TRUE,time_horizon = 1:3,learner = "learn_ranger"))
summary(x)
system.time(x <- run_rtmle(x,refit = TRUE,time_horizon = 2,folds = 5,learner = c("learn_glm","learn_glmnet","learn_ranger")))

######################################################################
### test-rtmle-machine-learning.R ends here
