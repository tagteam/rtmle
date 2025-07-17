### test-stochastic-dynamic-interventions.R ---
#----------------------------------------------------------------------
## Author: Alessandra Meddis
## Created: May 21st
## Version:
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
setwd("C:/Users/fkd749/Desktop/KVN2021/Amalie/rtmle-AM")
# source all the functions
for (f in list.files("R/",pattern = "R$",full.names = TRUE)){source(f)}


library(data.table)
#library(targets)
library(prodlim)
set.seed(17)
ld <- simulate_long_data(n = 11130,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))


protocol(x) <- list(name = "Always_A",
                    treatment_variables = "A",
                    intervention=intervene_dynamic,
                    intervention_type="dynamic")
prepare_data(x) <- list()
target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "Always_A")
x <- run_rtmle(x,refit = TRUE,time_horizon = 3,learner = "learn_glm")

summary(x)

set.seed(17)
ld <- simulate_long_data(n = 11130,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
y <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
y$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
add_baseline_data(y) <- ld$baseline_data[,start_followup_date:=0]
y <- long_to_wide(y,intervals = seq(0,2000,30.45*12))


protocol(y) <- list(name = "Always_A", intervention = data.frame("A" = factor("1",levels = c("0","1"))))
prepare_data(y) <- list()
target(y) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "Always_A")
y <- run_rtmle(y,refit = TRUE,time_horizon = 3,learner = "learn_glm")
summary(y)




protocol(y) <- list(name = "Always_A", treatment_variables="A",
                    intervention = intervene_stochastic,
                    intervention_type="stochastic")
prepare_data(y) <- list()
target(y) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "Always_A")
y <- run_rtmle(y,refit = TRUE,time_horizon = 3,learner = "learn_glm")
summary(y)



### by ltmle:
library(ltmle)
library(data.table)

ldata <- copy(x$prepared_data)
ldata[,id := NULL]
ldata[,rtmle_predicted_outcome := NULL]
ldata[,start_followup_date := NULL]
ldata[,A_3 := NULL]
ldata[,c("A_0","A_1","A_2") := lapply(.SD,function(a){1*(a == 1)}),.SDcols = c("A_0","A_1","A_2")]


res_l<-ltmle(data = ldata,
      Anodes = c("A_0","A_1","A_2"),
      Cnodes = c("Censored_1","Censored_2","Censored_3"),
      Ynodes = c("Y_1","Y_2","Y_3"),
      Lnodes = c("L_0","Dead_1","L_1","Dead_2","L_2"),
      gform = c(A_0 = "A_0 ~ sex + age + L_0",
      Censored_1 = "Censored_1 ~ sex + age + L_0 + A_0",A_1 = "A_1 ~ sex + age + L_0 + A_0",
      Censored_2 = "Censored_2 ~ sex + age + L_0 + A_0 + L_1 + A_1",A_2 = "A_2 ~ sex + age + L_0 + A_0 + L_1 + A_1",
      Censored_3 = "Censored_3 ~ sex + age + L_0 + A_0 + L_1 + A_1 + L_2 + A_2"),
      Qform = c(Y_1 = "Q.kplus1 ~ sex + age + L_0 + A_0",
                Y_2 = "Q.kplus1 ~ sex + age + L_0 + A_0 + L_1 + A_1",
                Y_3 = "Q.kplus1 ~ sex + age + L_0 + A_0 + L_1 + A_1 + L_2 + A_2"),
      survivalOutcome = TRUE,
      deterministic.Q.function = detQ,
      abar = rep(0,3),
      estimate.time = FALSE)

summary(res_l)
######################################################################
### test-stochastic-dynamic-interventions.R ends here





detQ <- function(data, current.node, nodes, called.from.estimate.g){
  death.index <- grep("Dead_",names(data))
  if(length(death.index)==0){
    ## message("No death/terminal event node found")
    return(NULL)
  }
  hist.death.index <- death.index[death.index < current.node]
  if(length(hist.death.index)==0)
    return(NULL)
  else{
    is.deterministic <- Reduce("+",lapply(data[,hist.death.index,drop=FALSE],
                                          function(dd){x=dd;x[is.na(dd)] <- 0;x}))>=1
    # should be unnecessary to exclude those who readily
    # have a missing value for death, but it does not hurt either
    is.deterministic[is.na(is.deterministic)] <- FALSE
    list(is.deterministic=is.deterministic, Q.value=0)
  }
}

