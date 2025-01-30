library(survival)
library(testthat)

test_that("sim_surv_data simulates data in the right way",{
  set.seed(857)
  beta <- matrix(rnorm(4,0,3), ncol = 2, nrow = 2)
  data <- sim_surv_data(N = 1000, beta = beta)

  survfit_death <- coxph(Surv(Time, Delta == 1) ~ I(L0 / 50) + A0, data = data)
  expect_true(confint(survfit_death)[1,1] <= beta[1,2] & beta[1,2] <= confint(survfit_death)[1,2])
  expect_true(confint(survfit_death)[2,1] <= beta[2,2] & beta[2,2] <= confint(survfit_death)[2,2])

  survfit_cens <- coxph(Surv(Time, Delta == 0) ~ I(L0 / 50) + A0, data = data)
  expect_true(confint(survfit_cens)[1,1] <= beta[1,1] & beta[1,1] <= confint(survfit_cens)[1,2])
  expect_true(confint(survfit_cens)[2,1] <= beta[2,1] & beta[2,1] <= confint(survfit_cens)[2,2])
})
