library(tidyverse)
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
#N <- 2000
#data_test <- sim_recur_event_data(N = N, beta = beta, eta = eta, nu = nu, at_risk = at_risk2, term_deltas = term_deltas)
## Creating a k and m variable
#data_test <- data_test %>% mutate(k = ave(ID, ID, FUN = seq_along) - 1) %>%
#  mutate(m = ave(ID, ID, FUN = seq_along))
#
#
#survfit_oper <- coxph(Surv(Time, Delta == 0) ~ L0 + k + A + L1, data = data_test[data_test$k == 0,])
#survfit_cens <- coxph(Surv(Time, Delta == 1) ~ L0 + k + A + L1, data = data_test[data_test$k == 0,])
#survfit_death <- coxph(Surv(Time, Delta == 2) ~ L0 + k + A + L1, data = data_test[data_test$k == 0,])
#survfit_cov <- coxph(Surv(Time, Delta == 3) ~ L0 + k + A + L1, data = data_test[data_test$k == 0,])
#confint(survfit_oper); beta0
#confint(survfit_cens); beta1
#confint(survfit_death); beta2
#confint(survfit_cov); beta3
