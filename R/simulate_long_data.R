### simulate_long_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 11 2024 (13:24) 
## Version: 
## Last-Updated: Sep 30 2024 (09:18) 
##           By: Thomas Alexander Gerds
##     Update #: 249
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#'
#' library(data.table)
#' ld <- simulate_long_data(n=100,number_epochs=20,register_format=TRUE)
#' map_grid(grid=seq(0,10*timevar_data,
#' ld[,table(table(id))]
#' ld[,num:=.N,by="id"]
#' ld[num==1]
#' ld
#' @export
simulate_long_data <- function(n,
                               number_epochs = 10,
                               baseline_rate,
                               beta,
                               register_format = FALSE,
                               int_dist=FALSE,
                               baseline_hazard_outcomes = 0.001) {
    beta_init <- list(A0_on_A = 0,
                      A0_on_Y = 0,
                      A0_on_D = 0,
                      A0_on_C = 0,
                      sum_A_on_Y = 0,
                      age_on_Y = 0,
                      age_on_D = 0,
                      age_on_C = 0,
                      age_on_A = 0,
                      age_on_L = 0,
                      sex_on_Y = 0,
                      sex_on_D = 0,
                      sex_on_C = 0,
                      sex_on_A = 0,
                      sex_on_L = 0,
                      L0_on_A = 0,
                      L0_on_Y = 0,
                      L0_on_D = 0,
                      L0_on_C = 0,
                      sum_L_on_Y = 0)
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
    # V = doctor visit time where treatment A is decided 
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
    if (int_dist){
      Baseline_Rate[["C"]] <- rep(0,nt)
    }
    
    # baseline variables
    pop <- data.table(
        id = 1:n,
        sex=rbinom(n,1,.4),
        age=runif(n,40,90),
        A = as.numeric(rep(NA, n)),
        sum_A = numeric(n),
        sum_L = numeric(n),
        time = numeric(n),
        event = rep("0",n)
    )
    # baseline treatment depends on baseline variables
    ## pop[, A_0:=rbinom(.N,1,lava::expit(0.35+0.1*L_0-0.3*sex-0.01*age))]
    if (int_dist)
      pop[, A_0:=1]
    else {
      pop[, A_0:=rbinom(.N,1,0.5)]
    }
  
    
    pop[, L_0:=rbinom(n,1,.17)]
    people_atrisk <- pop[,.(id,entrytime = time,age,sex,L_0,A_0,sum_L,A,sum_A)]

    if (int_dist){
      people_atrisk[,propensity_A := 1]
    } else {
      people_atrisk[,propensity_A := lava::expit(-3 + Beta$age_on_A*age+Beta$sex_on_A*sex+Beta$L0_on_A*L_0+Beta$A0_on_A*A_0)] ## does not depend on previous treatment; only baseline treatment
    }
    people_atrisk[,hazard_ratio_L := exp(Beta$age_on_L*age+Beta$sex_on_L*sex)]
    people_atrisk[,hazard_ratio_Y := exp(Beta$age_on_Y*age+Beta$sex_on_Y*sex+Beta$L0_on_Y*L_0+Beta$A0_on_Y*A_0)]
    people_atrisk[,hazard_ratio_D := exp(Beta$age_on_D*age+Beta$sex_on_D*sex+Beta$L0_on_D*L_0+Beta$A0_on_D*A_0)]
    people_atrisk[,hazard_ratio_C := exp(Beta$age_on_C*age+Beta$sex_on_C*sex+Beta$L0_on_C*L_0+Beta$A0_on_C*A_0)]
    # fup_info collects followup information has_terminal collects terminal information
    fup_info <- NULL
    has_terminal <- NULL
    # time loop
    j <- 1
    # doctor visit schedule (could be person specific but fixed for now)
    schedule <- seq(365,10*365,365)
    while (j < number_epochs && nrow(people_atrisk)>0){
        # calculate the time and type of the minimum of latent times to V,L,C,Y,D
        # matrix with latent times
        next_doctor_visit <- c(schedule,10*365)[1+prodlim::sindex(eval.times = people_atrisk$entrytime,jump.times = schedule)]
        ttt = do.call("cbind",list(
                                  next_doctor_visit,
                                  reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = Baseline_Rate[["L"]],hazardratio = people_atrisk$hazard_ratio_L,entrytime = people_atrisk$entrytime,decimals = 2),
                                  reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = Baseline_Rate[["C"]],hazardratio = people_atrisk$hazard_ratio_C,entrytime = people_atrisk$entrytime,decimals = 2),
                                  reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = Baseline_Rate[["Y"]],hazardratio = people_atrisk$hazard_ratio_Y,entrytime = people_atrisk$entrytime,decimals = 2),
                                  reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = Baseline_Rate[["D"]],hazardratio = people_atrisk$hazard_ratio_D,entrytime = people_atrisk$entrytime,decimals = 2))
                      )
        mins = Rfast::rowMins(ttt,value = FALSE)
        people_atrisk[,event := factor(mins,levels = 1:5,labels = c("V","L","C","Y","D"))]
        people_atrisk[,time := Rfast::rowMins(ttt,value = TRUE)]
        ## print(people_atrisk[id == 10,.(j = j,entrytime,time)])
        # censor at max_fup
        people_atrisk[time>max_fup,event := "C"]
        people_atrisk[time>max_fup,time := max_fup]
        is_terminal <- !(people_atrisk$event%in%c("V","L"))
        #------------------------------------------------------------------------------
        # collect terminal information
        #
        has_terminal <- rbind(has_terminal,
                              people_atrisk[is_terminal,.(id,terminal_time = time,terminal_event = event)])
        #------------------------------------------------------------------------------
        # restrict to people still at risk
        #
        people_atrisk = people_atrisk[!is_terminal]
        # draw treatment at the doctor visit times
        people_atrisk[event == "V",A := rbinom(.N,1,propensity_A)]

        # update propensity score
        # if (int_dist)
        #   people_atrisk[event == "V",propensity_A := 1]
        # else
        #   people_atrisk[event == "V",propensity_A := lava::expit(-3 + Beta$age_on_A*age+Beta$sex_on_A*sex+Beta$L0_on_A*L_0+Beta$A0_on_A*A)]
        people_atrisk[event == "V",sum_A := sum_A+A]
        # add to comorbidity index 
        people_atrisk[event == "L",sum_L := sum_L+1]
        # collect followup information
        fup_info <- rbind(fup_info,people_atrisk[,names(pop),with = FALSE],fill = TRUE)
        # -----------------------------------------------------------------------------
        # update for next epoch
        people_atrisk[,entrytime := time]
        people_atrisk[,hazard_ratio_L := hazard_ratio_L]
        people_atrisk[,hazard_ratio_Y := hazard_ratio_Y*exp(((Beta$sum_L_on_Y*sum_L)+(Beta$sum_A_on_Y*sum_A)))]
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
        bsl <- pop[time == 0,.(id,sex,age)]
        event_data <- pop[time == 0,.(id,terminal_time = terminal_time,terminal_event)]
        censored_data <- event_data[terminal_event == "C",.(id,date = terminal_time)]
        competing_data <- event_data[terminal_event == "D",.(id,date = terminal_time)]
        outcome_data <- event_data[terminal_event == "Y",.(id,date = terminal_time)]
        timevar_baseline <- pop[time == 0 & L_0 == 1,.(id,date = 0)]
        timevar_data <- pop[event == "L",.(id,date = time)]
        timevar_data <- rbind(timevar_baseline,timevar_data)
        # random baseline treatment
        # treatment_baseline <- pop[time == 0,.(id,date = 0, A = A_0)]
        # treatment_data <- pop[event == "V",.(id,date = time, A = A)]
        treatment_baseline <- pop[time == 0 & A_0 == 1,.(id,date = 0)]
        treatment_data <- pop[event == "V" & A == 1,.(id,date = time)] 
        treatment_data <- rbind(treatment_baseline,treatment_data)
        list(baseline_data = bsl,
             timevar_data = list(L = timevar_data,A = treatment_data),
             outcome_data = outcome_data,
             competing_data = competing_data,
             censored_data = censored_data)
    }else{
        pop[]
    }        
}



######################################################################
### simulate_long_data.R ends here
