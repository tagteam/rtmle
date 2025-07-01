library(testthat)
library(rtmle)
library(data.table)
library(prodlim)
test_that("run rtmle on simulated data",{
    set.seed(112)
    ld <- simulate_long_data(n = 1000,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(intervals = 3,
                    name_id = "id",
                    name_outcome = "Y",
                    name_competing = "Dead",
                    name_censoring = "Censored",
                    censored_label = "censored")
    x <- add_long_data(x,
                    outcome_data=ld$outcome_data,
                    censored_data=ld$censored_data,
                    competing_data=ld$competing_data,
                    timevar_data=ld$timevar_data)
    x <- add_baseline_data(x,data=ld$baseline_data)
    x <- long_to_wide(x,intervals = seq(0,2000,30.45*6))
    x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1)
    ## x <- prepare_data(x)
    x <- prepare_data(x) 
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "Always_A")
    x <- model_formula(x)
    suppressWarnings(x <- run_rtmle(x,verbose = FALSE))
    expect_equal(x$estimate$Main_analysis$Estimate,0.220619,tolerance = 0.001)
    expect_equal(x$estimate$Main_analysis$Standard_error,0.01502242,tolerance = 0.001)
})
test_that("obsolete way of run rtmle on simulated data",{
    set.seed(112)
    ld <- simulate_long_data(n = 1000,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(intervals = 2,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x <- add_long_data(x,outcome_data = ld[["outcome_data"]],censored_data = ld[["censored_data"]],competing_data = ld[["competing_data"]],timevar_data = ld[["timevar_data"]])
    x <- add_baseline_data(x,data = ld$baseline_data)
    x <- long_to_wide(x,intervals = seq(0,2000,30.45*6))
    x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1)
    x <- prepare_data(x) 
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "Always_A")
    x <- model_formula(x)
    suppressWarnings(x <- run_rtmle(x,verbose = FALSE))
    expect_equal(x$estimate$Main_analysis$Estimate,0.220619,tolerance = 0.001)
    expect_equal(x$estimate$Main_analysis$Standard_error,0.01502242,tolerance = 0.001)
})
test_that("run rtmle on without competing risks",{
    set.seed(112)
    ld <- simulate_long_data(n = 1000,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(intervals = 2,name_id = "id",name_outcome = "Y",name_competing = NULL,name_censoring = "Censored",censored_label = "censored")
    x$long_data <- ld[c("outcome_data","censored_data","timevar_data")]
    x$long_data$outcome_data <- rbind(x$long_data$outcome_data,ld$competing_data)
    setkey(x$long_data$outcome_data,id,date)
    x <- add_baseline_data(x,data=ld$baseline_data)
    x <- long_to_wide(x,intervals = seq(0,2000,30.45*6))
    x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1)
    x <- prepare_data(x)     
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "Always_A")
    x <- model_formula(x)
    suppressWarnings(x <- run_rtmle(x,verbose = FALSE))
    summary(x)
})



