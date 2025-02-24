### simulate_long_data.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds & Alessandra
## Created: Jul 11 2024 (13:24)
## Version:
## Last-Updated: Jan 30 2025 (11:53)
##           By: Alessandra
##     Update #: 269
#----------------------------------------------------------------------
##
### Commentary:  We want to simulate data with one OR two treatments variable (combination)
##               We want to generate data from a multi-logistic model
##               because we want the probability of getting "A", "A+B", "B", or none (ref)
### Change Log:
#----------------------------------------------------------------------
##
### Code:
#' Simulating longitudinal data for illustration purposes
#'
#' FIXME
#' @title Simulating longitudinal data
#' @param n Sample size
#' @param number_visits Number of doctor visit times covariates and treatment change
#' @param baseline_rate Vector of hazard rates
#' @param beta List of regression coefficients
#' @param register_format Logical. If \code{TRUE} the result is not in wide format but re-formatted as a list of register data.
#' @param interventional_distribution Logical. If \code{TRUE} the data generating mechanism produces uncensored data under the intervention
#' @param baseline_hazard_outcomes Baseline hazard function
#' @examples
#' simulate_long_data(10)
#' @export
simulate_long_data <- function(n,
                               number_visits = 10,
                               baseline_rate,
                               beta,
                               register_format = FALSE,
                               interventional_distribution=FALSE,
                               treatment_names="A",
                               baseline_hazard_outcomes = 0.001) {
  #hazard_ratio_C <- hazard_ratio_D <- hazard_ratio_L <- hazard_ratio_A <- hazard_ratio_Y <- event <- terminal_time <- terminal_event <- entrytime <- NULL
  #A_0 <- B_0<-L_0 <- id <- age <- sex <- sum_L <- A <-B<-sumb_B<- sum_A <- propensity_A <- NULL
  n_trt<-2
  ## We consider more treatments ( max 2) and these are the definition of the betas on the treatment(s)
  comb.grid <- expand.grid(treatment_names, treatment_names)
  treatment_betas<-as.list(rep(0,8*length(treatment_names) + nrow(comb.grid)))
  names(treatment_betas)<-c(paste0("",treatment_names,"0_on_A", sep=""),
                            paste0("",treatment_names,"0_on_Y", sep=""),
                            paste0("",treatment_names,"0_on_D", sep=""),
                            paste0("",treatment_names,"0_on_C", sep=""),
                            paste0("sum_",treatment_names,"on_Y", sep=""),
                            paste0("sum_",comb.grid[[1]],"_on_",comb.grid[[2]], sep=""),
                            paste0("age_on_",treatment_names,"", sep=""),
                            paste0("sex_on_",treatment_names,"", sep=""),
                            paste0("L0_on_",treatment_names,"", sep=""))

  # remaining betas:
  beta_init <- c(treatment_betas,list(age_on_Y = 0,
                                      age_on_D = 0,
                                      age_on_C = 0,
                                      age_on_L = 0,
                                      sex_on_Y = 0,
                                      sex_on_D = 0,
                                      sex_on_C = 0,
                                      sex_on_L = 0,
                                      L0_on_Y = 0,
                                      L0_on_D = 0,
                                      L0_on_C = 0,
                                      sum_L_on_Y = 0))


  if (missing(beta)) {
    Beta <- beta_init
  } else {
    Beta <- c(beta,beta_init)
    # remove the now obsolete init values
    Beta <- Beta[unique(names(Beta))]
  }
  # mimick a 10-year time frame
  time = seq(1,10*365,1)
  max_fup <- 10*365
  nt <- length(time)
  # V = doctor visit time where treatment A and/or B is decided
  # L = onset time of comorbidity (recurrent event)
  baseline_init = list(L = cumsum(rep(baseline_hazard_outcomes,nt)),
                       Y = cumsum(rep(baseline_hazard_outcomes,nt)),
                       D = cumsum(rep(baseline_hazard_outcomes,nt)),
                       C = cumsum(rep(baseline_hazard_outcomes,nt)))
  if (missing(baseline_rate)) {
    Baseline_Rate <- baseline_init
  } else {
    Baseline_Rate <- c(baseline_rate,baseline_init)
    # remove the now obsolete init values
    Baseline_Rate <- Baseline_Rate[unique(names(Baseline_Rate))]
  }
  if (interventional_distribution){
    Baseline_Rate[["C"]] <- rep(0,nt)
  }

  # baseline variables

  pop <- data.table(
    id = 1:n,
    sex=stats::rbinom(n,1,.4),
    age=stats::runif(n,40,90),
    sum_L = numeric(n),
    time = numeric(n),
    event = rep("0",n)
  )

  # add  variables about the treatment(s):
  vars_cumtrt<-paste0("sum_",treatment_names,"", sep="")
   pop[,(vars_cumtrt):=numeric(n)]
   pop[,(treatment_names):=as.numeric(rep(NA,n))]

  # baseline treatment depends on baseline variables
   vars_trt_baseline<-paste0("",treatment_names,"_0", sep="")
   pop[,(vars_trt_baseline):=as.numeric(rep(NA,n))]
  if (interventional_distribution)
    pop[,(vars_trt_baseline):=1]
  else {
    # we assume that both treatments are balanced, we could define this in a different way,
    # with some probability definition somewhere
    pop[, (vars_trt_baseline):=lapply(.SD, function(x){rbinom(.N,1,0.5)}) , .SDcols=vars_trt_baseline ]
  }

   # we will create a treatment variable that is the combination of the two
   # and this is a variable with 4 different values
   pop[, trt.comb:=rowSums(.SD), .SDcols=vars_trt_baseline]
   # do this in a more elegant way in the future
   if(length(treatment_names)>1){
   pop[A_0==0 & B_0==1, trt.comb:=3]
     n_trt=3
     }



  pop[, L_0:=stats::rbinom(n,1,.17)]
  pop$entrytime<-pop$time
  cols<-c("id","entrytime","age","sex","L_0","sum_L","trt.comb", vars_trt_baseline,vars_cumtrt,treatment_names)
  people_atrisk <- pop[, ..cols]

  ## we want a propensity score on the treatment combination and not on the single treatment
  if (interventional_distribution){
    people_atrisk[,propensity_trt := 1]
  } else {
    tmp.prop<-NULL
    for( r in 1:n_trt){
      ## if r=1--> we are calculating the predictors for A=1 and B=0
      ## if r=2--> we are calculating the predictors for A=1 and B=1
      # if we have only one treatment--> r=1 (only on time) and then it calculates it only for the one treatment
      tmp.trt.names=treatment_names[1:r]
      if(r==n_trt & r>1){ # we have 2 treatments (this has to change if we want more treatments)
      ## if r=3--> we are calculating the predictors for A=0 and B=1
      tmp.trt.names=treatment_names[-1]}
      ntrt<-length(tmp.trt.names)
      prop<-0
      ## linear predictor for the propensity score model:
     for(p in 1:ntrt){
        prop= prop + (-3 + rowSums( Beta[c(paste0("age_on_",tmp.trt.names[p],"",sep=""),paste0("sex_on_",tmp.trt.names[p],"",sep=""),
             paste0("L0_on_",tmp.trt.names[p],"",sep=""))] * people_atrisk[, .(age,sex,L_0) ] ))}

      tmp.prop<-cbind(tmp.prop,exp(prop))

    }

    #denomimator of the expit function for more categories
    den.prop<-1+apply(tmp.prop,1,sum)
    # calculate the propensity score for levels 1,2,3
    propensity_trt<-sapply(1:n_trt, FUN=function(s){tmp.prop[,s]/(den.prop) } )
    #-- add the one for trt=0
    propensity_trt<-cbind(1-rowSums(propensity_trt),propensity_trt)
    prop_names<-paste0("propensity_comb",0:max(pop$trt.comb))
    #we have to add this o people_at_risk so tp use if later:
    people_atrisk[,(prop_names):=lapply( (0:n_trt) +1 ,FUN=function(x){propensity_trt[,x]})]
    #people_atrisk[,propensity_trt := lava::expit(-3 + Beta$age_on_A*age+Beta$sex_on_A*sex+Beta$L0_on_A*L_0)]
  }
  ###################################
  ## maybe there is a better way to do this with more treatments, directly inside the data frame
  tot_trt_Y<-rowSums(Beta[paste0("",treatment_names,"0_on_Y")] * people_atrisk[, ..vars_trt_baseline])
  tot_trt_D<-rowSums(Beta[paste0("",treatment_names,"0_on_D")] * people_atrisk[, ..vars_trt_baseline])
  tot_trt_C<-rowSums(Beta[paste0("",treatment_names,"0_on_C")] * people_atrisk[, ..vars_trt_baseline])

  people_atrisk[,hazard_ratio_L := exp(Beta$age_on_L*age+Beta$sex_on_L*sex)]
  people_atrisk[,hazard_ratio_Y := exp(Beta$age_on_Y*age+Beta$sex_on_Y*sex+Beta$L0_on_Y*L_0+tot_trt_Y)]
  people_atrisk[,hazard_ratio_D := exp(Beta$age_on_D*age+Beta$sex_on_D*sex+Beta$L0_on_D*L_0+tot_trt_D)]
  people_atrisk[,hazard_ratio_C := exp(Beta$age_on_C*age+Beta$sex_on_C*sex+Beta$L0_on_C*L_0+tot_trt_C)]
  # fup_info collects followup information has_terminal collects terminal information
  fup_info <- NULL
  has_terminal <- NULL
  # time loop
  j <- 1
  # doctor visit schedule (could be person specific but fixed for now)
  schedule <- seq(365,10*365,365)
  while (j < number_visits && nrow(people_atrisk)>0){
    # calculate the time and type of the minimum of latent times to V,L,C,Y,D
    # matrix with latent times
    next_doctor_visit <- c(schedule,10*365)[1+prodlim::sindex(eval.times = people_atrisk$entrytime,jump.times = schedule)]
    ttt = do.call("cbind",list(
      next_doctor_visit,
      reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = Baseline_Rate[["L"]],
                 hazardratio = people_atrisk$hazard_ratio_L,entrytime = people_atrisk$entrytime,decimals = 2),
      reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = Baseline_Rate[["C"]],
                 hazardratio = people_atrisk$hazard_ratio_C,entrytime = people_atrisk$entrytime,decimals = 2),
      reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = Baseline_Rate[["Y"]],
                 hazardratio = people_atrisk$hazard_ratio_Y,entrytime = people_atrisk$entrytime,decimals = 2),
      reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = Baseline_Rate[["D"]],
                 hazardratio = people_atrisk$hazard_ratio_D,entrytime = people_atrisk$entrytime,decimals = 2))
    )
    mins = Rfast::rowMins(ttt,value = FALSE)
    people_atrisk[,event := factor(mins,levels = 1:5,labels = c("V","L","C","Y","D"))]
    people_atrisk[,time := Rfast::rowMins(ttt,value = TRUE)]
    ## print(people_atrisk[id == 10,data.table::data.table(j = j,entrytime,time)])
    # censor at max_fup
    people_atrisk[time>max_fup,event := "C"]
    people_atrisk[time>max_fup,time := max_fup]
    is_terminal <- !(people_atrisk$event%in%c("V","L"))
    #------------------------------------------------------------------------------
    # collect terminal information
    #
    has_terminal <- rbind(has_terminal,
                          people_atrisk[is_terminal,data.table::data.table(id,terminal_time = time,terminal_event = event)])
    #------------------------------------------------------------------------------
    # restrict to people still at risk
    #
    people_atrisk = people_atrisk[!is_terminal]
    # draw treatment at the doctor visit times
    # for more treatments then we have the propensity score for the different combinations 0,1,2,3
    # 0: A=0, B=0
    # 1: A=1,B=0
    # 2: A=1,B=1
    # 3: A=0, B=1
    people_atrisk[event == "V", trt.comb:=sample(0:n_trt, size=1, prob=.SD),.SDcols=prop_names, by=id ]

    ### Update the treatment values A, B at this time, absed on the trt.combination
    if(n_trt>2){
    people_atrisk[event == "V", (treatment_names) := .(fcase(
      trt.comb==0, 0,
      trt.comb==1, 1,
      trt.comb==2, 1,
      trt.comb==3, 0),
      fcase(
        trt.comb==0, 0,
        trt.comb==1, 0,
        trt.comb==2, 1,
        trt.comb==3, 01)) ]

      }
    else{
      people_atrisk[event == "V", (treatment_names) :=trt.comb]
    }
    # update the cumulative treatment that are part of the propensity score after the first time:
    people_atrisk[event == "V", (vars_cumtrt) := lapply(1:length(treatment_names), FUN=function(x){get(vars_cumtrt[x]) + get(treatment_names[x])})]


    # update propensity score for the next time:
    if (interventional_distribution)
      people_atrisk[event == "V",propensity_A := 1] ## we have to change this cause it is not the correct intervention
    else{
    tmp.prop<-NULL
    for( r in 1:n_trt){

      if(r==n_trt & r>1){
        tmp.trt.names=treatment_names[-1]}
      else{
        tmp.trt.names=treatment_names[1:r]
      }
      ntrt<-length(tmp.trt.names)
      prop<-0
      for(p in 1:ntrt){
        prop= prop + (-3 + rowSums( Beta[c(paste0("age_on_",tmp.trt.names[p],"",sep=""),paste0("sex_on_",tmp.trt.names[p],"",sep=""),
                                           paste0("L0_on_",tmp.trt.names[p],"",sep=""),
                                           paste0("sum_",treatment_names,"_on_",tmp.trt.names[p],sep=""))] * people_atrisk[event=="V", .SD, .SDcols=c("age","sex","L_0",treatment_names)] ))}

      tmp.prop<-cbind(tmp.prop,exp(prop))
    }

    den.prop<-1+apply(tmp.prop,1,sum)
    propensity_trt<-sapply(1:n_trt, FUN=function(s){tmp.prop[,s]/(den.prop) } )
    propensity_trt<-cbind(1-rowSums(propensity_trt),propensity_trt)
    people_atrisk[event=="V",(prop_names):=lapply( (0:n_trt) +1 ,FUN=function(x){propensity_trt[,x]})]
    }

    # add to comorbidity index
    people_atrisk[event == "L",sum_L := sum_L+1]
    # collect followup information
    fup_info <- rbind(fup_info,people_atrisk[,names(pop),with = FALSE],fill = TRUE)
    # -----------------------------------------------------------------------------
    # update for next epoch
    people_atrisk[,entrytime := time]
    people_atrisk[,hazard_ratio_L := hazard_ratio_L]
    # cumulative treatment for the outcome
    tot_cumtrt_Y<-rowSums(Beta[paste0("",vars_cumtrt,"on_Y")] * people_atrisk[, ..vars_cumtrt])
    people_atrisk[,hazard_ratio_Y := hazard_ratio_Y*exp((Beta$sum_L_on_Y*sum_L) + tot_cumtrt_Y) ]
    people_atrisk[,hazard_ratio_D := hazard_ratio_D]
    people_atrisk[,hazard_ratio_C := hazard_ratio_C]
    j = j+1
  }
  # combine terminal information with followup
  pop <- has_terminal[pop,on = "id"]
  fup_info <- has_terminal[fup_info,on = "id"]
  pop <- rbind(pop,fup_info)
  # if missing terminal event and reached max_fup then censor
  pop[is.na(terminal_event) & time == max_fup,`:=`(terminal_time = max_fup,terminal_event = "C")]
  # if still any missing terminal event draw one
  if (length(needs_terminal <- unique(pop[is.na(terminal_event)][["id"]]))>0){
    p <- people_atrisk[data.table(id = needs_terminal),on = "id"]
    ttt = do.call("cbind",list(
      reventtime(n = nrow(p),breaks = time,cumhazard = Baseline_Rate[["C"]],hazardratio = p$hazard_ratio_C,entrytime = people_atrisk$entrytime,decimals = 2),
      reventtime(n = nrow(p),breaks = time,cumhazard = Baseline_Rate[["Y"]],hazardratio = p$hazard_ratio_Y,entrytime = people_atrisk$entrytime,decimals = 2),
      reventtime(n = nrow(p),breaks = time,cumhazard = Baseline_Rate[["D"]],hazardratio = p$hazard_ratio_D,entrytime = people_atrisk$entrytime,decimals = 2))
    )
    mins = Rfast::rowMins(ttt,value = FALSE)
    p[,terminal_event := factor(mins,levels = 1:3,labels = c("C","Y","D"))]
    p[,terminal_time := Rfast::rowMins(ttt,value = TRUE)]
    pop[is.na(terminal_event),`:=`(terminal_time = 10*365,terminal_event = "C")]
  }
  setkey(pop,id,time,event)

  # clean up for those who directly have a terminal event
  pop[is.na(event),event := terminal_event]
  pop[is.na(time),time := terminal_time]
  pop[,terminal_event := factor(terminal_event,levels = c("Y","D","C"))]

  if (register_format){
    bsl <- pop[time == 0,data.table::data.table(id,sex,age)]
    event_data <- pop[time == 0,data.table::data.table(id,terminal_time = terminal_time,terminal_event)]
    censored_data <- event_data[terminal_event == "C",data.table::data.table(id,date = terminal_time)]
    competing_data <- event_data[terminal_event == "D",data.table::data.table(id,date = terminal_time)]
    outcome_data <- event_data[terminal_event == "Y",data.table::data.table(id,date = terminal_time)]
    timevar_baseline <- pop[time == 0 & L_0 == 1,data.table::data.table(id,date = 0)]
    timevar_data <- pop[event == "L",data.table::data.table(id,date = time)]
    timevar_data <- rbind(timevar_baseline,timevar_data)
    ## treatment baseline: (at least one is 1 and if only 1 treatment-->trt.comb=0)
    # but we need one data for data set:
    treatment_list<-as.list(rep(0,length(treatment_names)))
    names(treatment_list)<-treatment_names
    for(t in 1:length(treatment_names)){
      treatment_baseline <- pop[time == 0 & get(vars_trt_baseline[t])==1,data.table::data.table(id,date = 0)]
      treatment_time <- pop[event == "V"  & get(treatment_names[t])==1,data.table::data.table(id,date = time)]
      treatment_list[[t]] <- rbind(treatment_baseline,treatment_time)
    }


    ## we have to define time_vardata depending on the different treatments


    list(baseline_data = bsl,
         timevar_data = c(list(L = timevar_data),treatment_list),
         outcome_data = outcome_data,
         competing_data = competing_data,
         censored_data = censored_data)
  }else{
    pop[]
  }
}



######################################################################
### simulate_long_data.R ends here
