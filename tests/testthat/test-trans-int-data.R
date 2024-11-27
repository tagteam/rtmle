library(testthat)
test_that("trans-int-data transforms data correctly",{
  set.seed(856)
  data1 <- sim_op_data(50)
  data2 <- trans_int_data(data1)
  expect_true(nrow(data1) == nrow(data2))
})
