library(survival)
library(testthat)

test_that("sim_op_data simulates data in the right way",{
  set.seed(858)
  # Generate data
  beta <- matrix(rnorm(16,0,1), ncol = 4, nrow = 4)
  data_test <- sim_op_data(500, beta = beta)

  # Transform data into tstart tstop format
  data_int <- trans_int_data(data_test)

  # Defining at risk indicators
  data_int[, at_risk_oper := as.numeric(A == 0)]
  data_int[, at_risk_cov := as.numeric(L == 0)]


  # Fit models
  survfit_cens <- coxph(Surv(tstart, tstop, Delta == 0) ~ I(L0 / 50) + A0 + L + A, data = data_int, cluster = ID)
  survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ I(L0 / 50) + A0 + L + A, data = data_int, cluster = ID)
  survfit_oper <- coxph(Surv(tstart, tstop, Delta == 2) ~ I(L0 / 50) + A0 + L, data = data_int[at_risk_oper == 1], cluster = ID)
  survfit_cov <- coxph(Surv(tstart, tstop, Delta == 3) ~ I(L0 / 50) + A0 + A, data = data_int[at_risk_cov == 1], cluster = ID)

  # Compare confidence intervals and true values
  expect_true(confint(survfit_oper, level = 0.99)[1,1] <= beta[1,3] & beta[1,3] <= confint(survfit_oper, level = 0.99)[1,2])
  expect_true(confint(survfit_oper, level = 0.99)[2,1] <= beta[2,3] & beta[2,3] <= confint(survfit_oper, level = 0.99)[2,2])
  expect_true(confint(survfit_oper, level = 0.99)[3,1] <= beta[3,3] & beta[3,3] <= confint(survfit_oper, level = 0.99)[3,2])
  expect_true(confint(survfit_death, level = 0.99)[1,1] <= beta[1,2] & beta[1,2] <= confint(survfit_death, level = 0.99)[1,2])
  expect_true(confint(survfit_death, level = 0.99)[2,1] <= beta[2,2] & beta[2,2] <= confint(survfit_death, level = 0.99)[2,2])
  expect_true(confint(survfit_death, level = 0.99)[3,1] <= beta[3,2] & beta[3,2] <= confint(survfit_death, level = 0.99)[3,2])
  expect_true(confint(survfit_death, level = 0.99)[4,1] <= beta[4,2] & beta[4,2] <= confint(survfit_death, level = 0.99)[4,2])
  expect_true(confint(survfit_cens, level = 0.99)[1,1] <= beta[1,1] & beta[1,1] <= confint(survfit_cens, level = 0.99)[1,2])
  expect_true(confint(survfit_cens, level = 0.99)[2,1] <= beta[2,1] & beta[2,1] <= confint(survfit_cens, level = 0.99)[2,2])
  expect_true(confint(survfit_cens, level = 0.99)[3,1] <= beta[3,1] & beta[3,1] <= confint(survfit_cens, level = 0.99)[3,2])
  expect_true(confint(survfit_cens, level = 0.99)[4,1] <= beta[4,1] & beta[4,1] <= confint(survfit_cens, level = 0.99)[4,2])
  expect_true(confint(survfit_cov, level = 0.99)[1,1] <= beta[1,4] & beta[1,4] <= confint(survfit_cov, level = 0.99)[1,2])
  expect_true(confint(survfit_cov, level = 0.99)[2,1] <= beta[2,4] & beta[2,4] <= confint(survfit_cov, level = 0.99)[2,2])
  expect_true(confint(survfit_cov, level = 0.99)[3,1] <= beta[4,4] & beta[4,4] <= confint(survfit_cov, level = 0.99)[3,2])
})
