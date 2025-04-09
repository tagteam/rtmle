### test-subset-rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr  8 2025 (14:47) 
## Version: 
## Last-Updated: Apr  9 2025 (11:30) 
##           By: Thomas Alexander Gerds
##     Update #: 14
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
test_that("run rtmle on a subset",{
    tau <- 2
    set.seed(17)
    ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
    add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
    x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
    protocol(x) <- list(name = "Always_A",intervention = data.frame("A" = factor("1",levels = c("0","1"))),verbose = FALSE)
    protocol(x) <- list(name = "Never_A",intervention = data.frame("A" = factor("0",levels = c("0","1"))),verbose = FALSE)
    prepare_data(x) <- list()
    target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = c("Always_A","Never_A"))
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
    x1$long_data <- ld1[c("outcome_data","censored_data","competing_data","timevar_data")]
    add_baseline_data(x1) <- ld1$baseline_data[,start_followup_date:=0]
    x1 <- long_to_wide(x1,intervals = seq(0,2000,30.45*12))
    protocol(x1) <- list(name = "Always_A",intervention = data.frame("A" = factor("1",levels = c("0","1"))),verbose = FALSE)
    protocol(x1) <- list(name = "Never_A",intervention = data.frame("A" = factor("0",levels = c("0","1"))),verbose = FALSE)
    prepare_data(x1) <- list()
    target(x1) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = c("Always_A","Never_A"))
    x1 <- run_rtmle(x1,learner = "learn_glmnet",time_horizon = 2,verbose = FALSE)
    x79 <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 2,subsets = list(S = list(label = "S79",id = 1:79)),verbose = FALSE)
    ## rbind(x1$estimate[[1]][[1]],x79$estimate$S79[[1]][[1]])
    expect_equal(x1$estimate,x79$estimate$S79)
})

test_that("stratified analyses",{
    tau <- 2
    set.seed(17)
    ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
    add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
    x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
    protocol(x) <- list(name = "Always_A",intervention = data.frame("A" = factor("1",levels = c("0","1"))),verbose = FALSE)
    protocol(x) <- list(name = "Never_A",intervention = data.frame("A" = factor("0",levels = c("0","1"))),verbose = FALSE)
    prepare_data(x) <- list()
    target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = c("Always_A","Never_A"))
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau,verbose = FALSE)
    # stratified analyses
    sex_strata <- list(list(label="Sex",append = TRUE,variable="Sex",value="Female",id=x$prepared_data[sex==0,id]),
                                list(label="Sex",append = TRUE,variable="Sex",value="Male",id=x$prepared_data[sex==1,id]))
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = tau,verbose=FALSE,subsets=sex_strata)
})

######################################################################
### test-subset-rtmle.R ends here
