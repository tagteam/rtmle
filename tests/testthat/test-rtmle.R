library(testthat)
library(rtmle)
library(data.table)
library(prodlim)
set.seed(112)
ld <- simulate_long_data(n = 1000,number_epochs = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
x <- rtmle_init(intervals = 2,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
x <- long_to_wide(x,intervals = seq(0,2000,30.45*6))
protocol(x) <- list(name = "Always_A",
                    treatment_variables = "A",
                    intervention = 1)
prepare_data(x) <- list(reset = TRUE)
target(x) <- list(name = "Outcome_risk",
                  strategy = "additive",
                  estimator = "tmle",
                  protocols = "Always_A")
suppressWarnings(x <- run_rtmle(x))
expect_equal(x$estimate$Outcome_risk$Always_A$Estimate,0.220619,tolerance = 0.001)
expect_equal(x$estimate$Outcome_risk$Always_A$Standard_error,0.01502242,tolerance = 0.001)



