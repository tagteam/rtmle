### testing.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 25 2024 (09:49) 
## Version: 
## Last-Updated: Nov 16 2024 (16:51) 
##           By: Thomas Alexander Gerds
##     Update #: 5
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
library(data.table)
library(prodlim)
library(targets)
setwd("~/research/SoftWare/rtmle/")

parse_learners(c(list("learn_glmnet",
                      "glm"=list(learn_variables="A")),
                 list("learn_ranger"),
                 list(list(learner_fun="learn_ranger",num.trees=5))))

#tar_source("R/")
set.seed(112)
ld <- simulate_long_data(n = 10000,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6))
set.seed(18124)
ld0 <- simulate_long_data(n = 10000,number_visits = 20,beta = list(A0_on_Y = 0,A0_on_A = 6))
ld[,table(terminal_event)]
ld[,table(A_0)]
ld[,table(table(id))]

plot(prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,reverse = TRUE,data = ld[!duplicated(id)]),cause = "Y",xlim = c(0,365.25*3),plot.main = "Outcome risk, protective effect of A on outcome",confint = 0L)

par(mfrow = c(2,2))
plot(prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,data = ld[!duplicated(id)]),cause = "Y",xlim = c(0,365.25*3),plot.main = "Outcome risk, protective effect of A on outcome",confint = 0L)
plot(prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,data = ld[!duplicated(id)]),cause = "D",xlim = c(0,365.25*3),plot.main = "Risk of death without outcome, protective effect of A on outcome",confint = 0L)
plot(prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,data = ld0[!duplicated(id)]),cause = "Y",xlim = c(0,365.25*3),plot.main = "Outcome risk, no effect of A",confint = 0L)
plot(prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,data = ld0[!duplicated(id)]),cause = "D",xlim = c(0,365.25*3),plot.main = "Risk of death without outcome, no effect of A",confint = 0L)



######################################################################
### testing.R ends here
