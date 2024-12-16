###### RTMLE  # CTP data in longitudinal form
# Jeg forsøger at følge manualen i rtmle
time_var_data <- readRDS('./sandbox/CTP/time_var_data.rds')
baseline_data <- readRDS('./sandbox/CTP/baseline_data.rds')
outcome_dt <- readRDS("./sandbox/CTP/outcome_dt.rds")
treatment_dt <- readRDS("./sandbox/CTP/treatment_dt.rds")
setnames(outcome_dt,"pnr","ID")
setnames(treatment_dt,"pnr","ID")

library(rtmle)
library(data.table)

x <- rtmle_init(intervals=5,name_id='ID',name_time='period',name_outcome='outcome',name_competing='compete',name_censoring='censor',treatment_levels=c('0','1'),censored_levels=c('1','0'),censored_label="1")
x$data <- list(baseline_data = baseline_data[,.(ID,Age)],
               outcome_data = outcome_dt,
               timevar_data = treatment_dt)
prepare_data(x) <- list()
protocol(x) <- list(name = "treat",treatment_variables = "treatment",intervention = 1)
protocol(x) <- list(name = "placebo",treatment_variables = "treatment",intervention = 0)
target(x) <- list(name = "Risk",strategy = "additive",estimator = "tmle",protocols = c("treat","placebo"))
x <- run_rtmle(x,time_horizon = 4)
summary(x)


################################################
################################################

# CTP data i wide format, klar til analyse
# Her viser jeg den gamle prepareTMLE selv om den ikke er aktuel mere,
# men så kan man se hvad der er hvad.
# Her kan jeg ikke gennemskue hvordan jeg skal forklare rtmle-objektet hvad der er hvad
outcome_dt <- readRDS("./sandbox/CTP/outcome_dt.rds")
treatment_dt <- readRDS("./sandbox/CTP/treatment_dt.rds")
## LTMLE
age <- baseline_data[,.(ID,Age)]
setnames(age,"ID","pnr")
time_horizon <- 5
x <- prepare_Ltmle(
  outcome_data = outcome_dt,
  regimen_data = treatment_dt,
  baseline_data = age, #need baseline covariates
  timevar_data = NULL, #need timevarying covariates
  time_horizon = 1:time_horizon,
  censored_label = "1",
  name_outcome = "outcome",
  name_regimen = "treatment",
  name_censoring = "censor",
  name_competing_risk = "compete",
  abar = rep(1,time_horizon),
  SL.library = "glm",
  verbose = TRUE,
  gbounds = c(0, 1),
  deterministic.g.function = MaintainControl ## if you are not on treatment, then you      #will not be on the treatment at a later data
)
fit_treat <- do.call("Ltmle",x)
x$abar <- rep(0,time_horizon)
fit_ref <- do.call("Ltmle",x)
