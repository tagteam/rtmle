### test-stochastic-dynamic-interventions.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  1 2024 (10:37) 
## Version: 
## Last-Updated: Nov  1 2024 (10:54) 
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
library(rtmle)
library(data.table)
library(targets)
library(prodlim)
set.seed(17)
ld <- simulate_long_data(n = 11130,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
my_break <- function(data,time){
    browser(skipCalls=1L)
    if (time == 0)
        return(1)
    else
        data[,{.SD},.SDcols = paste0("A_",0:(time-1))]
}
protocol(x) <- list(name = "Always_A",treatment_variables = "A",intervention = "my_break")
prepare_data(x) <- list()
target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "Always_A")
system.time(x <- run_rtmle(x,refit = TRUE,time_horizon = 2,learner = "learn_glm"))

######################################################################
### test-stochastic-dynamic-interventions.R ends here
