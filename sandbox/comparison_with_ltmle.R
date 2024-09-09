### comparison_with_ltmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 25 2024 (09:50) 
## Version: 
## Last-Updated: Aug 22 2024 (16:07) 
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
# compare with ltmle
library(rtmle)
library(data.table)
library(targets)
library(prodlim)
tar_source("~/research/Methods/TMLE_for_breakfast/Ltmle/R/")
source("~/research/Epidemi/Reddie/LEADER/functions/run_ltmle.R")
source("~/research/Epidemi/Reddie/LEADER/functions/summary.runLtmle.R")
set.seed(112)
ld <- simulate_long_data(n = 100,number_epochs = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)

x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
add_long_data(x) <- ld
prepare_data(x) <- list(treatment_variables = "A",reset = TRUE,intervals = seq(0,2000,30.45*6))
protocol(x) <- list(name = "Always_A",treatment_variables = "A",intervention = 1)
target(x) <- list(name = "Outcome_risk",strategy = "additive", estimator = "tmle", estimands = 3,protocol = "Always_A")
system.time(x <- run_rtmle(x))
# Ltmle
vn <- names(x$prepared_data)
w_treatment <- x$prepared_data[,c("id",grep("A_",vn,value = TRUE)),with = FALSE]
w_outcome <- x$prepared_data[,c("id",grep("Y_|Censored_|Dead_",vn,value = TRUE)),with = FALSE]
w_timevar <- x$prepared_data[,c("id",grep("L_",vn,value = TRUE)),with = FALSE]
system.time(tfit <- run_ltmle(name_outcome="Y",time_horizon=3,reduce = FALSE,regimen_data=list("A" = w_treatment),outcome_data=list("Y" = w_outcome),baseline_data=x$prepared_data[,.(id,sex,age)],timevar_data=w_timevar,SL.library="glm",censor_others = FALSE,gbounds=c(0,1),abar = rep(1,3),name_id = "id",verbose=FALSE,gcomp = FALSE))
## system.time(tfit1 <- run_ltmle(name_outcome="Y",time_horizon=3,reduce = FALSE,regimen_data=list("A" = w_treatment),outcome_data=list("Y" = w_outcome),baseline_data=x$prepared_data[,.(id,sex,age)],timevar_data=w_timevar,SL.library="glm",censor_others = FALSE,gbounds=c(0,1),abar = list(rep(1,3),rep(0,3)),name_id = "id",verbose=FALSE,gcomp = TRUE))
summary(tfit)
source("~/research/SoftWare/rtmle/R/summary.rtmle.R")
summary(x)
## summary(tfit1)
all.equal(as.numeric(tfit$A$Ltmle_fit$estimate),x$targets$Outcome_risk$estimate)
all.equal(c(tfit$A$Ltmle_fit$IC),x$targets$Outcome_risk$IC)



######################################################################
### comparison_with_ltmle.R ends here
