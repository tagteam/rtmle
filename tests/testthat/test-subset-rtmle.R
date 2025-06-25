### test-subset-rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr  8 2025 (14:47) 
## Version: 
## Last-Updated: Jun 18 2025 (07:08) 
##           By: Thomas Alexander Gerds
##     Update #: 31
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
    tau <- 2
    set.seed(17)
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
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau,verbose = FALSE)
    # subset to first 79
    ld1 <- lapply(ld,function(d){
        subset(d,d$id <= 79)
    })
    names(ld1) <- names(ld)
    ld1$timevar_data <- lapply(ld$timevar_data,function(d){
        subset(d,d$id <= 79)
    })
    names(ld1$timevar_data) <- names(ld$timevar_data)
    x1 <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x1 <- add_long_data(x1,
                    outcome_data=ld1$outcome_data,
                    censored_data=ld1$censored_data,
                    competing_data=ld1$competing_data,
                    timevar_data=ld1$timevar_data)
    x1 <- add_baseline_data(x1,data = ld1$baseline_data)
    x1 <- long_to_wide(x1,intervals = seq(0,2000,30.45*12))
    x1 <- protocol(x1,name = "Always_A",intervention = data.frame("A" = factor("1",levels = c("0","1"))),verbose = FALSE)
    x1 <- protocol(x1,name = "Never_A",intervention = data.frame("A" = factor("0",levels = c("0","1"))),verbose = FALSE)
    x1 <- prepare_data(x1)
    x1 <- target(x1,name = "Outcome_risk",estimator = "tmle",protocols = c("Always_A","Never_A"))
    x1 <- model_formula(x1)
    x1 <- run_rtmle(x1,learner = "learn_glmnet",time_horizon = 2,verbose = FALSE)
    x79 <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 2,subsets = list(S = list(label = "S79",id = 1:79)),verbose = FALSE)
    ## rbind(x1$estimate[[1]][[1]],x79$estimate$S79[[1]][[1]])
    data.table::setattr(x79$estimate$S79,"IC",NULL)
    expect_equal(x1$estimate$Main_analysis,x79$estimate$S79)
})

test_that("stratified analyses",{
    tau <- 2
    set.seed(17)
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
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau,verbose = FALSE)
    # stratified analyses
    sex_strata <- list(list(label="Sex",append = TRUE,variable="Sex",level="Female",id=x$prepared_data[sex==0,id]),
                                list(label="Sex",append = TRUE,variable="Sex",level="Male",id=x$prepared_data[sex==1,id]))
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = tau,verbose=FALSE,subsets=sex_strata)
})

######################################################################
### test-subset-rtmle.R ends here
