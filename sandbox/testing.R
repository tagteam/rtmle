### testing.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 25 2024 (09:49) 
## Version: 
## Last-Updated: maj 28 2026 (13:06) 
##           By: Thomas Alexander Gerds
##     Update #: 16
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

library(rtmle)
library(data.table)

data(simulated_cohort)
ld <- register_format(simulated_cohort)
x <- rtmle_init(time_grid = seq(0, 20, 4),name_id = "id",name_outcome = "stroke",name_competing = "death",name_censoring = "dropout",censored_label = "censored")
x <- add_long_data(x,outcome_data = ld$timevar_data$stroke[!duplicated(id)],censored_data = ld$timevar_data$dropout,competing_data = ld$timevar_data$death,timevar_data = ld$timevar_data[c("bleeding", "changeSBP", "A", "B")])
x <- add_baseline_data(x, data = ld$baseline_data)
x <- long_to_wide(x, start_followup_date = 0)
x <- prepare_rtmle_data(x)
x <- protocol(x,name = "Use_A_not_B_in_the_end",intervention = data.frame(time = x$intervention_nodes,A = factor(c(NA,NA,NA,rep("1",2)),levels = c("0", "1")),B = factor(c(NA,NA,NA,"0","0"),levels = c("0", "1"))),verbose = FALSE)
x <- protocol(x,name = "Always_A_not_B_in_the_beginning",intervention = data.frame(time = x$intervention_nodes,A = factor(rep("1",5),levels = c("0", "1")),B = factor(c("0","0",NA,NA,NA),levels = c("0", "1"))),verbose = FALSE)
x <- target(x,name = "Stroke_risk",protocols = c("Use_A_not_B_in_the_end", "Always_A_not_B_in_the_beginning"),estimator = "tmle")
x <- model_formula(x, verbose = FALSE)
x <- run_rtmle(x, learner = "learn_glm", time_horizon = 1:2, verbose = FALSE)
summary(x)


parse_learners(c(list("learn_glmnet",
                      "glm"=list(learn_variables="A")),
                 list("learn_ranger"),
                 list(list(fun="learn_ranger",num.trees=5))))


set.seed(112)
ld <- simulate_long_data(n = 1000,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),baseline_hazard_outcomes = 0.001,register_format = TRUE)
x <- rtmle_init(time_grid = seq(0,2000,30.45*6),name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
x <- add_long_data(x,outcome_data=ld$outcome_data,censored_data=ld$censored_data,competing_data=ld$competing_data,timevar_data=ld$timevar_data)
x <- add_baseline_data(x,data=ld$baseline_data)
x <- long_to_wide(x,start_followup_date = 0)
x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1)
x <- prepare_rtmle_data(x)
x$prepared_data[Y_2 == 1,Y_2 := 0]
# test if constant variable is removed
x$prepared_data[,L_1 := 1]
x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "Always_A")
x <- model_formula(x)
x <- run_rtmle(x,verbose = FALSE,learner = "learn_glmnet")
x$models
summary(x)
plot_model_parameters_all(x, Model_outcomes=c("Y"))


dt <- rbindlist(list(ld$timevar_data$A[,event := "A"],
                     ld$timevar_data$L[,event := "L"],
                     ld$outcome_data[,event := "Y"],
                     ld$competing_data[,event := "D"],
                     ld$censored_data[,event := "C"]))
setkey(dt,id,date)



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
