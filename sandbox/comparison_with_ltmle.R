### comparison_with_ltmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 25 2024 (09:50) 
## Version: 
## Last-Updated: Jan 10 2025 (14:48) 
##           By: Thomas Alexander Gerds
##     Update #: 55
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
## source("~/research/Epidemi/Reddie/LEADER/functions/run_ltmle.R")
## source("~/research/Epidemi/Reddie/LEADER/functions/summary.runLtmle.R")
set.seed(17)
tau <- 3
ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
protocol(x) <- list(name = "Always_A",treatment_variables = "A",intervention = 1)
prepare_data(x) <- list()
target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "Always_A")
x <- run_rtmle(x,learner = "learn_glm",time_horizon = 1:tau)
## summary(x)
# Ltmle
vn <- names(x$prepared_data)
w_treatment <- x$prepared_data[,c("id",grep("A_",vn,value = TRUE)),with = FALSE]
w_outcome <- x$prepared_data[,c("id",grep("Y_|Censored_|Dead_",vn,value = TRUE)),with = FALSE]
w_timevar <- x$prepared_data[,c("id",grep("L_",vn,value = TRUE)),with = FALSE]
tfit <- run_ltmle(name_outcome="Y",time_horizon=1:tau,reduce = FALSE,regimen_data=list("A" = w_treatment),outcome_data=list("Y" = w_outcome),baseline_data=x$prepared_data[,.(id,sex,age)],timevar_data=w_timevar,SL.library="glm",censor_others = FALSE,gbounds=c(0,1),abar = rep(1,tau),name_id = "id",verbose=FALSE,gcomp = FALSE)
tfit1 <- run_ltmle(stratify = TRUE,
                   name_outcome="Y",
                   time_horizon=1:tau,
                   reduce = FALSE,
                   regimen_data=list("A" = w_treatment),
                   outcome_data=list("Y" = w_outcome),
                   baseline_data=x$prepared_data[,.(id,sex,age)],
                   timevar_data=w_timevar,
                   SL.library="glm",
                   censor_others = FALSE,
                   gbounds=c(0,1),
                   abar = rep(1,tau),
                   name_id = "id",
                   verbose=FALSE,
                   gcomp = FALSE)
summary(tfit)
## summary(x)
all.equal(as.numeric(tfit$A$Ltmle_fit$estimate),x$estimate$Outcome_risk[["Always_A"]]$Estimate)
all.equal(c(tfit$A$Ltmle_fit$IC),unlist(x$IC$Outcome_risk$Always_A,use.names = FALSE))
## tfit$A$Ltmle_fit$fit$Q
## tfit$A$Ltmle_fit$fit$Qstar
## b <- tfit$A$Ltmle_fit$formulas$Qform
## x$models
## a <- sapply(x$models[["Always_A"]][["outcome"]],function(x)x$formula)


library(ltmle)
ldata <- copy(x$prepared_data)
ldata[,id := NULL]
ldata[,rtmle_predicted_outcome := NULL]
ldata[,start_followup_date := NULL]
y <- ltmle(data = ldata,
           Anodes = c("A_0","A_1","A_2"),
           Cnodes = c("Censored_1","Censored_2","Censored_3"),
           Ynodes = c("Y_1","Y_2","Y_3"),
           Lnodes = c("L_0","Dead_1","L_1","Dead_2","L_2"),
           survivalOutcome = TRUE,
           deterministic.Q.function = detQ,
           abar = rep(1,3),
           estimate.time = FALSE,
           )
summary(y)
summary(tfit$A$Ltmle_fit)


if (FALSE){
    # formulas
    a <- tfit$A$Ltmle_fit$formulas
    A <- y$formulas
    all.equal(a,A)
    b <- sapply(x$models[["Always_A"]],function(x)x$formula)
    cbind(c(a$gform,a$Qform)[names(b)],b)
    # fitted g-models
    for (g in c("A_0","Censored_1","A_1","Censored_2","A_2","Censored_3")){
        print(g)
        a <- tfit$A$Ltmle_fit$fit$g[[g]]
        b <- data.frame(x$models[["Always_A"]][[g]]$fit,check.names = FALSE)
        print(all.equal(a,b))
    }
    # fitted Q-models
    for (q in c("Y_1","Y_2","Y_3")){
        print(q)
        a <- tfit$A$Ltmle_fit$fit$Q[[q]]
        b <- data.frame(x$models[["Always_A"]][[q]]$fit,check.names = FALSE)
        print(all.equal(a,b))
    }
}


######################################################################
### comparison_with_ltmle.R ends here
