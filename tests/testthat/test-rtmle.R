library(testthat)
library(rtmle)
library(data.table)
library(prodlim)
library(targets)
set.seed(112)
ld <- simulate_long_data(n = 10000,number_epochs = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
x <- rtmle_init(intervals = 8,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
add_long_data(x) <- ld
prepare_data(x) <- list(treatment_variables = "A",
                        reset = TRUE,
                        intervals = seq(0,2000,30.45*6))
protocol(x) <- list(name = "Always_A",
                    treatment_variables = "A",
                    intervention = 1)
target(x) <- list(name = "Outcome_risk", strategy = "additive",
                  estimator = "tmle", estimands = 3, protocol = "Always_A")
system.time(x <- run_rtmle(x))



