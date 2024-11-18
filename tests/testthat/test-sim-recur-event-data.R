library(dplyr)
library(survival)
set.seed(856)

# Survival and competing at risk function
# You are only at risk until some event happens
at_risk1 <- function(x, k, m) as.numeric(k == 0)

# Operation setting at risk function
at_risk2 <- function(x, k, m) {
  # If you have not died yet or been censored yet, you are at risk for dying or being censored
  if(x == 1 | x == 2) return(1)
  # You are only at risk for an operation if you have not had an operation yet
  else if(x == 0) return(as.numeric(k == 0 | (k == 1 & m == 1)))
  # You are only at risk for a change in the covariate process if you have not had a change yet
  else return(as.numeric(m == 0))
}


# Testing whether estimated effects correspond to true effects

# Generating data
N <- 2000
data_test <- sim_recur_event_data(N = 2000)
# Creating a k and m variable
data_test <- data_test %>% mutate(k = ave(ID, ID, FUN = seq_along) - 1) #%>%
  #mutate(m = ave(ID, ID, FUN = seq_along))


survfit_oper <- coxph(Surv(Time, Delta == 0) ~ L0 + k + A + L1, data = data_test[data_test$k == 0,])
survfit_cens <- coxph(Surv(Time, Delta == 1) ~ L0 + k + A + L1, data = data_test[data_test$k == 0,])
survfit_death <- coxph(Surv(Time, Delta == 2) ~ L0 + k + A + L1, data = data_test[data_test$k == 0,])
survfit_cov <- coxph(Surv(Time, Delta == 3) ~ L0 + k + A + L1, data = data_test[data_test$k == 0,])
confint(survfit_oper)[1,1] <= 0 & 0 <= confint(survfit_oper)[1,2]
confint(survfit_oper)[3,1] <= 0 & 0 <= confint(survfit_oper)[3,2]
confint(survfit_oper)[4,1] <= 0 & 0 <= confint(survfit_oper)[4,2]
confint(survfit_cens)[1,1] <= 0 & 0 <= confint(survfit_cens)[1,2]
confint(survfit_cens)[3,1] <= 0 & 0 <= confint(survfit_cens)[3,2]
confint(survfit_cens)[4,1] <= 0 & 0 <= confint(survfit_cens)[4,2]
confint(survfit_death)[1,1] <= 0 & 0 <= confint(survfit_death)[1,2]
confint(survfit_death)[3,1] <= 0 & 0 <= confint(survfit_death)[3,2]
confint(survfit_death)[4,1] <= 0 & 0 <= confint(survfit_death)[4,2]
confint(survfit_cov)[1,1] <= 0 & 0 <= confint(survfit_cov)[1,2]
confint(survfit_cov)[3,1] <= 0 & 0 <= confint(survfit_cov)[3,2]
confint(survfit_cov)[4,1] <= 0 & 0 <= confint(survfit_cov)[4,2]
