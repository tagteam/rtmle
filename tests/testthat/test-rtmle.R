library(testthat)
library(rtmle)
library(data.table)
library(prodlim)
test_that("run rtmle on simulated data",{
    set.seed(112)
    ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(time_grid = seq(0,1500,30.45*6),name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x <- add_long_data(x,
                       outcome_data=ld$outcome_data,
                       censored_data=ld$censored_data,
                       competing_data=ld$competing_data,
                       timevar_data=ld$timevar_data)
    x <- add_baseline_data(x,data=ld$baseline_data)
    x <- long_to_wide(x,start_followup_date=0)
    x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1)
    x <- prepare_rtmle_data(x) 
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "Always_A")
    x <- model_formula(x)
    suppressWarnings(x <- run_rtmle(x,time_horizon = 2,refit = TRUE,verbose=0L))
    expect_equal(x$estimate$Main_analysis$Estimate,0.1497072,tolerance = 0.001)
    expect_equal(x$estimate$Main_analysis$Standard_error,0.06217111,tolerance = 0.001)
})

test_that("censoring models are reused across protocols",{
    set.seed(42)
    n <- 120
    W <- rnorm(n)
    A_0 <- factor(rbinom(n,1,plogis(W)),levels = c(0,1))
    C_1 <- factor(ifelse(rbinom(n,1,plogis(-0.5 + 0.8 * W + 0.7 * (A_0 == "1"))) == 1,
                         "censored",
                         "uncensored"),
                  levels = c("uncensored","censored"))
    Y_1 <- rbinom(n,1,plogis(-1 + W + 0.5 * (A_0 == "1")))
    Y_1[C_1 == "censored"] <- NA
    x <- rtmle_init(time_grid = 0:1,
                    name_id = "id",
                    name_outcome = "Y",
                    name_competing = NULL,
                    name_censoring = "C",
                    censored_label = "censored")
    x$prepared_data <- data.table(id = 1:n,W = W,A_0 = A_0,C_1 = C_1,Y_1 = Y_1)
    x$names$name_baseline_covariates <- "W"
    x$names$name_time_covariates <- "A"
    x <- protocol(x,name = "Always_A",
                  intervention = data.frame(time = x$intervention_nodes,
                                            A = factor("1",levels = c("0","1"))),
                  verbose = FALSE)
    x <- protocol(x,name = "Never_A",
                  intervention = data.frame(time = x$intervention_nodes,
                                            A = factor("0",levels = c("0","1"))),
                  verbose = FALSE)
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = c("Always_A","Never_A"))
    x <- model_formula(x,verbose = FALSE)
    fit_calls <- 0
    reuse_calls <- 0
    counting_glm <- function(character_formula,
                             data,
                             intervened_data,
                             save_fitted_objects = FALSE,
                             reuse_fit = NULL,
                             ...){
        if (is.null(reuse_fit)){
            fit_calls <<- fit_calls + 1
        }else{
            reuse_calls <<- reuse_calls + 1
        }
        learn_glm(character_formula = character_formula,
                  data = data,
                  intervened_data = intervened_data,
                  save_fitted_objects = save_fitted_objects,
                  reuse_fit = reuse_fit,
                  ...)
    }
    x <- run_rtmle(x,
                   learner = list(name = "counting_glm",fun = counting_glm),
                   time_horizon = 1,
                   verbose = FALSE)
    expect_equal(fit_calls,5)
    expect_equal(reuse_calls,1)
    expect_true(is.list(x$models$time_0$censoring$C_1$fit))
})


test_that("run rtmle without covariates",{
    set.seed(112)
    ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(time_grid = seq(0,1500,30.45*6),name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    #  add treatment variables but no timevarying covariates
    x <- add_long_data(x,outcome_data=ld$outcome_data,censored_data=ld$censored_data,competing_data=ld$competing_data,timevar_data = ld$timevar_data["A"])
    #  need id variable in the baseline data
    x <- add_baseline_data(x,data=ld$baseline_data[,.(id)])
    x <- long_to_wide(x,start_followup_date=0)
    x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1)
    x <- prepare_rtmle_data(x) 
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "Always_A")
    x <- model_formula(x)
    suppressWarnings(x <- run_rtmle(x,time_horizon = 2,refit = TRUE,verbose = FALSE))
    expect_equal(x$estimate$Main_analysis$Estimate,0.1388889,tolerance = 0.001)
    expect_equal(x$estimate$Main_analysis$Standard_error,0.05795775,tolerance = 0.001)
})

test_that("run rtmle without competing risks",{
    set.seed(112)
    ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(time_grid = seq(0,1500,30.45*6),name_id = "id",name_outcome = "Y",name_competing = NULL,name_censoring = "Censored",censored_label = "censored")
    x$long_data <- ld[c("outcome_data","censored_data","timevar_data")]
    x$long_data$outcome_data <- rbind(x$long_data$outcome_data,ld$competing_data)
    setkey(x$long_data$outcome_data,id,date)
    x <- add_baseline_data(x,data=ld$baseline_data)
    x <- long_to_wide(x,start_followup_date = 0)
    x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1)
    x <- prepare_rtmle_data(x,verbose = FALSE)     
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "Always_A")
    x <- model_formula(x)
    suppressWarnings(x <- run_rtmle(x,verbose = FALSE))
    expect_output(print(summary(x)))
})

test_that("rtmle can use g-formula",{
    set.seed(112)
    ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(time_grid = seq(0,1500,30.45*6),name_id = "id",name_outcome = "Y",name_competing = NULL,name_censoring = "Censored",censored_label = "censored")
    x$long_data <- ld[c("outcome_data","censored_data","timevar_data")]
    x$long_data$outcome_data <- rbind(x$long_data$outcome_data,ld$competing_data)
    setkey(x$long_data$outcome_data,id,date)
    x <- add_baseline_data(x,data=ld$baseline_data)
    x <- long_to_wide(x,start_followup_date = 0)
    x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1)
    x <- prepare_rtmle_data(x,verbose = FALSE)     
    x <- target(x,name = "Outcome_risk",protocols = "Always_A")
    x <- model_formula(x)
    x <- run_rtmle(x,estimator = "tmle",verbose = FALSE)
    x$protocols$Always_A$cumulative_intervention_probs <- pmax(x$protocols$Always_A$cumulative_intervention_probs,0.5)
    z <- run_rtmle(x,refit = FALSE,estimator = "tmle",verbose = FALSE)
    y <- run_rtmle(x,estimator = "g-formula",verbose = FALSE)
    expect_output(print(summary(x)))    
})
