library(testthat)
library(rtmle)
library(data.table)
library(prodlim)
library(targets)
set.seed(126) 
base_hazard_outcome <- 0.0001
ld <- simulate_long_data(n = 10000,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 2.6),register_format = TRUE, baseline_hazard_outcomes = base_hazard_outcome)
time_horizon <- 9
res <- rep(NA, time_horizon)
for (time in 1:time_horizon){
  print(paste0("time: ", time))
  x <- rtmle_init(intervals = time,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
  add_long_data(x) <- ld
  protocol(x) <- list(name = "Always_A",
                      treatment_variables = "A",
                      intervention = 1)
  
  prepare_data(x) <- list(treatment_variables = "A",
                          reset = TRUE,
                          intervals = seq(0,10*365,365))
  target(x) <- list(name = "Outcome_risk", strategy = "additive",
                    estimator = "tmle", estimands = 3, protocols = "Always_A")
  capture.output(x <- run_rtmle(x))
  res[time] <- x$estimate$Outcome_risk$Always_A
}

ld_int <- simulate_long_data(200000,number_visits=20,beta=list(A_on_Y = -.2,A0_on_Y = -0.3),int_dist=TRUE,baseline_hazard_outcomes = base_hazard_outcome)
ld_int <- ld_int[, c("id", "terminal_time", "terminal_event")]
ld_int <- unique(ld_int)  

res_df <- data.table(
  time = 1:time_horizon * 365,
  risk_ltmle = res,
  risk_int = ld_int[, sapply(1:time_horizon, function(time)
    mean(terminal_time <= time * 365 & terminal_event == "Y"))]
)

library(ggplot2)
ggplot(res_df, aes(x = time)) +
  geom_line(aes(y = risk_ltmle, color = "LTMLE")) +
  geom_line(aes(y = risk_int, color = "True value (continuous time)")) +
  scale_color_manual(values = c("LTMLE" = "red", "True value (continuous time)" = "blue")) +
  labs(title = "Cumulative incidence function",
       x = "Time (days)",
       y = "Cumulative incidence") +
  theme_bw() +
  theme(legend.title = element_blank())
