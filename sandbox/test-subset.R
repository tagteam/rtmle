### test-subset.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: sep 20 2026 (08:08) 
## Version: 
## Last-Updated: sep 20 2026 (08:11) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
data(simulated_cohort)
ld <- register_format(simulated_cohort)
x <- rtmle_init(time_grid = seq(0,20,4),name_id = "id",name_outcome = "stroke",name_competing = "death",name_censoring = "dropout",censored_label = "censored")
x <- add_long_data(x,outcome_data=ld$timevar_data$stroke[!duplicated(id)],censored_data=ld$timevar_data$dropout,competing_data=ld$timevar_data$death,timevar_data=ld$timevar_data[c("bleeding","changeSBP","A","B")])
x <- add_baseline_data(x,data=ld$baseline_data)
x <- discretize_data(x,start_followup_date=0)
x <- prepare_rtmle_data(x)
x <- regime(x,name = "Always_A",intervention = data.frame(time=x$intervention_nodes,"A" = factor("1",levels = c("0","1"))))
x <- regime(x,name = "Never_A",intervention = data.frame(time=x$intervention_nodes,"A" = factor("0",levels = c("0","1"))))
x <- target(x,name = "Outcome_risk",estimator = "tmle",regimes = c("Always_A","Never_A"))
x <- model_formula(x)
x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 3)
# stratified analyses
x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 3,
               verbose=FALSE,
               subsets=list(list(label="Sex",variable="Sex",
                                 level="Female",id=x$prepared_data[sex==0,id]),
                            list(label="Sex",variable="Sex",
                                 level="Male",id=x$prepared_data[sex==1,id])))
x$estimate
######################################################################
### test-subset.R ends here
