### comparison_with_ltmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 25 2024 (09:50) 
## Version: 
## Last-Updated: Mar 17 2025 (14:12) 
##           By: Thomas Alexander Gerds
##     Update #: 76
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
breakfast_code <- "~/research/Methods/TMLE_for_breakfast/Ltmle/R/"
if (file.exists(breakfast_code)){
    tar_source(breakfast_code)
} else{
    stop(paste0("Please change the variable 'breakfast_code' to point into the",
                " right folder on your computer.\nYou can always download the folder here:",
                "https://github.com/tagteam/TMLE_for_breakfast/tree/main/Ltmle"))
}

# ------------------------------------------------------------------------------------------
# Intervening on single treatment variable
# ------------------------------------------------------------------------------------------
set.seed(17)
tau <- 3
ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
protocol(x) <- list(name = "Always_A",intervention = data.frame("A" = factor(1,levels = c(0,1))))
prepare_data(x) <- list()
target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "Always_A")
x <- run_rtmle(x,learner = "learn_glm",time_horizon = 1:tau)

# Ltmle
vn <- names(x$prepared_data)
w_treatment <- x$prepared_data[,c("id",grep("A_",vn,value = TRUE)),with = FALSE]
w_treatment[,c("A_0","A_1","A_2","A_3") := lapply(.SD,function(a){1*(a == 1)}),.SDcols = c("A_0","A_1","A_2","A_3")]
w_outcome <- x$prepared_data[,c("id",grep("Y_|Censored_|Dead_",vn,value = TRUE)),with = FALSE]
w_timevar <- x$prepared_data[,c("id",grep("L_",vn,value = TRUE)),with = FALSE]
tfit <- run_ltmle(name_outcome="Y",time_horizon=1:tau,reduce = FALSE,regimen_data=list("A" = w_treatment),outcome_data=list("Y" = w_outcome),baseline_data=x$prepared_data[,.(id,sex,age)],timevar_data=w_timevar,SL.library="glm",censor_others = FALSE,gbounds=c(0,1),abar = rep(1,tau),name_id = "id",verbose=FALSE,gcomp = FALSE)
tfit1 <- run_ltmle(stratify = TRUE,name_outcome="Y",time_horizon=1:tau,reduce = FALSE,regimen_data=list("A" = w_treatment),outcome_data=list("Y" = w_outcome),baseline_data=x$prepared_data[,.(id,sex,age)],timevar_data=w_timevar,SL.library="glm",censor_others = FALSE,gbounds=c(0,1),abar = rep(1,tau),name_id = "id",verbose=FALSE,gcomp = FALSE)
summary(tfit)
summary(x)
all.equal(as.numeric(tfit$A$Ltmle_fit$estimate),x$estimate$Outcome_risk[["Always_A"]]$Estimate)
all.equal(c(tfit$A$Ltmle_fit$IC),unlist(x$IC$Outcome_risk$Always_A,use.names = FALSE))
tfit$A$Ltmle_fit$fit$Q
tfit$A$Ltmle_fit$fit$Qstar
b <- tfit$A$Ltmle_fit$formulas$Qform
x$models
formulas <- sapply(x$models[["Always_A"]],function(x)x$formula)
# ------------------------------------------------------------------------------------------
# Intervening on more than 1 treatment variable
# ------------------------------------------------------------------------------------------

set.seed(17)
tau <- 2
ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
# generate a second treatment B at random
ld$timevar_data$B <- ld$timevar_data$A[,.(id = sort(sample(1:91,size = 75,replace = FALSE)),date = date)]
x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
protocol(x) <- list(name = "Always_A_Never_B",intervention = data.frame("A" = factor(1,levels = c(0,1)),"B" = factor(0,levels = c(0,1))))
prepare_data(x) <- list()
# modifying the prepared data 
x$prepared_data[,B_0 := rep(0,.N)]
x$names$name_constant_variables <- c("B_0",x$names$name_constant_variables)
target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "Always_A_Never_B")
x <- run_rtmle(x,refit = TRUE,learner = "learn_glm",time_horizon = tau)
# Ltmle
vn <- names(x$prepared_data)
w_treatment <- x$prepared_data[,c("id",grep("A_",vn,value = TRUE)),with = FALSE]
w_treatment[,c("A_0","A_1","A_2") := lapply(.SD,function(a){1*(a == 1)}),.SDcols = c("A_0","A_1","A_2")]
w_other <- x$prepared_data[,c("id",grep("B_",vn,value = TRUE)),with = FALSE]
w_other[,c("B_0","B_1","B_2") := lapply(.SD,function(a){1*(a == 1)}),.SDcols = c("B_0","B_1","B_2")]
w_outcome <- x$prepared_data[,c("id",grep("Y_|Censored_|Dead_",vn,value = TRUE)),with = FALSE]
w_timevar <- x$prepared_data[,c("id",grep("L_",vn,value = TRUE)),with = FALSE]
ab_treatment <- w_treatment[w_other,on = "id"]
## setcolorder(ab_treatment,c("id","A_0","B_0","A_1","B_1","A_2","B_2"))
setcolorder(ab_treatment,c("id","A_0","B_0","A_1","B_1"))
tfit <- run_ltmle(name_outcome="Y",time_horizon=tau,reduce = FALSE,regimen_data=list("A" = ab_treatment),outcome_data=list("Y" = w_outcome),baseline_data=x$prepared_data[,.(id,sex,age)],timevar_data=w_timevar,SL.library="glm",censor_others = TRUE,gbounds=c(0,1),name_id = "id",verbose=FALSE,gcomp = FALSE)

## unlist(tfit$A$Ltmle_fit$formulas)
## sapply(x$models[[1]],function(x)x$formula)
uu <- head(tfit$A$Ltmle_fit$cum.g.used)
Lcum <- tfit$A$Ltmle_fit$cum.g[uu]
Rcum <- x$cumulative_intervention_probs[[1]][,-2][uu]
all.equal(Lcum,Rcum)

colnames(x$cumulative_intervention_probs[[1]])
head(x$intervention_probs[[1]])[,-c(1,3)]

all.equal(tfit$A$Ltmle_fit$fit$g$A_1,
          x$models$Always_A_Never_B$A_1$fit)


summary(tfit)
summary(x)

c(LTMLE = as.numeric(tfit$A$Ltmle_fit$estimate),
  RTMLE = x$estimate$Outcome_risk[["Always_A_Never_B"]]$Estimate)
all.equal(as.numeric(tfit$A$Ltmle_fit$estimate),x$estimate$Outcome_risk[["Always_A_Never_B"]]$Estimate)
all.equal(c(tfit$A$Ltmle_fit$IC),unlist(x$IC$Outcome_risk$Always_A_Never_B,use.names = FALSE))

# FixedTimeTMLE can output the intervention_match matrix
## a <- 1*readRDS(file = "~/tmp/a.rds")
## b <- x$intervention_match[[1]]
## head(a)
## head(b)

# ------------------------------------------------------------------------------------------
# Comparison with CRAN-ltmle, that is without run_ltmle
# ------------------------------------------------------------------------------------------
if (FALSE){
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
