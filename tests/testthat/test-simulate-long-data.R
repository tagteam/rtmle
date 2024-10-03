### test-simulate-long-data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 24 2024 (13:05) 
## Version: 
## Last-Updated: Oct  2 2024 (15:59) 
##           By: Thomas Alexander Gerds
##     Update #: 13
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
library(rtmle)
library(prodlim)
library(data.table)
# no effect of A
set.seed(18124)
ld0 <- simulate_long_data(n = 100000,number_epochs = 20,beta = list(A0_on_A = 6),register_format = FALSE)
# effect of baseline A but no additional effect of A when used after baseline
set.seed(18124)
ldA0 <- simulate_long_data(n = 100000,number_epochs = 20,beta = list(A0_on_Y = -0.3,A0_on_A = 6))
# no effect of baseline A but effect of A when used after baseline 
set.seed(18124)
ldA <- simulate_long_data(n = 100000,number_epochs = 20,beta = list(sum_A_on_Y = -0.5,A0_on_A = 6))
# effect of baseline A and additional effect of A when used after baseline
set.seed(18124)
ldA0A <- simulate_long_data(n = 100000,number_epochs = 20,beta = list(sum_A_on_Y = -0.5,A0_on_Y = -0.3,A0_on_A = 6))
# effect of starting A 
fit0 <- prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,data = ld0[!duplicated(id)])
fitA0 <- prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,data = ldA0[!duplicated(id)])
fitA <- prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,data = ldA[!duplicated(id)])
fitA0A <- prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,data = ldA0A[!duplicated(id)])
# effect on Y
# FIXME: reactivate these tests
## expect_equal(diff(as.numeric(unlist(predict(fit0,times = 800,cause = "Y",newdata = data.frame(A_0 = c(0,1)))))),0,tolerance = 0.01)
## expect_equal(diff(as.numeric(unlist(predict(fitA0,times = 800,cause = "Y",newdata = data.frame(A_0 = c(0,1)))))),-0.0842,tolerance = 0.001)
## expect_equal(diff(as.numeric(unlist(predict(fitA,times = 800,cause = "Y",newdata = data.frame(A_0 = c(0,1)))))),-0.0621,tolerance = 0.001)
## expect_equal(diff(as.numeric(unlist(predict(fitA0A,times = 800,cause = "Y",newdata = data.frame(A_0 = c(0,1)))))),-0.135,tolerance = 0.001)

######################################################################
### test-simulate-long-data.R ends here
