library(survival)
library(testthat)

#test_that("sim_op_data simulates data in the right way",{
#  set.seed(857)
#  # Generate data
#  beta <- matrix(rnorm(16,0,1), ncol = 4, nrow = 4)
#  data_test <- sim_op_data(300, beta = beta)
#
#  # Transform data into tstart tstop format
#  data_int <- trans_int_data(data_test)
#
#  # Fit models
#  survfit_oper <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + k + A, data = data_int, cluster = ID)
#  survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + k + A + L, data = data_int, cluster = ID)
#  survfit_cens <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + k + A + L, data = data_int, cluster = ID)
#  survfit_cov <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0 + k + A, data = data_int, cluster = ID)
#
#  # Compare confidence intervals and true values
#  expect_true(confint(survfit_oper)[1,1] <= beta[1,1] & beta[1,1] <= confint(survfit_oper)[1,2])
#  #expect_true(confint(survfit_oper)[2,1] <= beta[2,1] + beta[4,1] & beta[2,1] + beta[4,1] <= confint(survfit_oper)[2,2])
#  expect_true(confint(survfit_oper)[3,1] <= beta[3,1] & beta[3,1] <= confint(survfit_oper)[3,2])
#  expect_true(confint(survfit_death)[1,1] <= beta[1,2] & beta[1,2] <= confint(survfit_death)[1,2])
#  expect_true(confint(survfit_death)[2,1] <= beta[2,2] & beta[2,2] <= confint(survfit_death)[2,2])
#  expect_true(confint(survfit_death)[3,1] <= beta[3,2] & beta[3,2] <= confint(survfit_death)[3,2])
#  expect_true(confint(survfit_death)[4,1] <= beta[4,2] & beta[4,2] <= confint(survfit_death)[4,2])
#  expect_true(confint(survfit_cens)[1,1] <= beta[1,3] & beta[1,3] <= confint(survfit_cens)[1,2])
#  expect_true(confint(survfit_cens)[2,1] <= beta[2,3] & beta[2,3] <= confint(survfit_cens)[2,2])
#  expect_true(confint(survfit_cens)[3,1] <= beta[3,3] & beta[3,3] <= confint(survfit_cens)[3,2])
#  expect_true(confint(survfit_cens)[4,1] <= beta[4,3] & beta[4,3] <= confint(survfit_cens)[4,2])
#  expect_true(confint(survfit_cov)[1,1] <= beta[1,4] & beta[1,4] <= confint(survfit_cov)[1,2])
#  expect_true(confint(survfit_cov)[2,1] <= beta[2,4] & beta[2,4] <= confint(survfit_cov)[2,2])
#  expect_true(confint(survfit_cov)[3,1] <= beta[3,4] & beta[3,4] <= confint(survfit_cov)[3,2])
#})
