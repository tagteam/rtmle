library(survival)
library(testthat)

test_that("sim_cr_data simulates data in the right way",{
  set.seed(856)
  beta <- matrix(rnorm(6,0,3), ncol = 3, nrow = 2)
  data <- sim_cr_data(N = 300, beta = beta)

  survfit_proc1 <- coxph(Surv(Time, Delta == 0) ~ L0 + A0, data = data)
  expect_true(confint(survfit_proc1)[1,1] <= beta[1,1] & beta[1,1] <= confint(survfit_proc1)[1,2])
  expect_true(confint(survfit_proc1)[2,1] <= beta[2,1] & beta[2,1] <= confint(survfit_proc1)[2,2])

  survfit_proc2 <- coxph(Surv(Time, Delta == 1) ~ L0 + A0, data = data)
  expect_true(confint(survfit_proc2)[1,1] <= beta[1,2] & beta[1,2] <= confint(survfit_proc2)[1,2])
  expect_true(confint(survfit_proc2)[2,1] <= beta[2,2] & beta[2,2] <= confint(survfit_proc2)[2,2])

  survfit_cens <- coxph(Surv(Time, Delta == 2) ~ L0 + A0, data = data)
  expect_true(confint(survfit_cens)[1,1] <= beta[1,3] & beta[1,3] <= confint(survfit_cens)[1,2])
  expect_true(confint(survfit_cens)[2,1] <= beta[2,3] & beta[2,3] <= confint(survfit_cens)[2,2])
})
