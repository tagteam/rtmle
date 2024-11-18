### test-rtmle-versus-ltmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov 16 2024 (17:04) 
## Version: 
## Last-Updated: Nov 18 2024 (11:29) 
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
library(testthat)
library(ltmle)
library(rtmle)
library(data.table)

testthat("compare with ltmle: not longitudinal",
{
    rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
    # Single time point Example
    n <- 1000
    W <- rnorm(n)
    A <- rexpit(-1 + 2 * W)
    Y <- rexpit(W + A)
    data <- data.frame(W, A, Y)
    result1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1)
    summary(result1)
    x <- rtmle_init(intervals = 1,name_id = "id",name_outcome = "Y",name_competing = NULL,name_censoring = NULL)
    rdata <- cbind(id = 1:n,data)
    setDT(rdata)
    setnames(rdata,c("id","W","A_0","Y_1"))
    x$prepared_data <- rdata
    x$names$name_baseline_covariates <- "W"
    x$names$name_time_covariates <- "A"
    protocol(x) <- list(name = "A",treatment_variables = "A",intervention = 1)
    target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "A")
    x <- run_rtmle(x,refit = TRUE)
    expect_equal(result1$estimates[[1]],x$estimate$Outcome_risk$A$Estimate)
})

testthat("compare with ltmle",{
    set.seed(17)
    ld <- simulate_long_data(n = 11130,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
    add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
    x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
    protocol(x) <- list(name = "Always_A",treatment_variables = "A",intervention = 1)
    prepare_data(x) <- list()
    target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "Always_A")
    x <- run_rtmle(x,learner = "learn_glm",time_horizon = 3)
    summary(x)
    ltmle()
})




######################################################################
### test-rtmle-versus-ltmle.R ends here
