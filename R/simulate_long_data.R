### simulate_long_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 11 2024 (13:24) 
## Version: 
## Last-Updated: Mar 25 2025 (12:35) 
##           By: Thomas Alexander Gerds
##     Update #: 274
#----------------------------------------------------------------------
## 
### Commentary: 
## 
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
                               baseline_hazard_outcomes = 0.001) {
    hazard_ratio_C <- hazard_ratio_D <- hazard_ratio_L <- hazard_ratio_A <- hazard_ratio_Y <- event <- terminal_time <- terminal_event <- entrytime <- NULL
    A_0 <- L_0 <- id <- age <- sex <- sum_L <- A <- sum_A <- propensity_A <- NULL
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
    if (interventional_distribution){
      Baseline_Rate[["C"]] <- rep(0,nt)
    }
    
    # baseline variables
    pop <- data.table(
        id = 1:n,
        sex=stats::rbinom(n,1,.4),
        age=stats::runif(n,40,90),
        A = as.numeric(rep(NA, n)),
        sum_A = numeric(n),
        sum_L = numeric(n),
        time = numeric(n),
        event = rep("0",n)
    )
    # baseline treatment depends on baseline variables
    ## pop[, A_0:=stats::rbinom(.N,1,lava::expit(0.35+0.1*L_0-0.3*sex-0.01*age))]
    if (interventional_distribution)
      pop[, A_0:=1]
    else {
      pop[, A_0:=stats::rbinom(.N,1,0.5)]
    }
  
    
    pop[, L_0:=stats::rbinom(n,1,.17)]
    people_atrisk <- pop[,data.table::data.table(id,entrytime = time,age,sex,L_0,A_0,sum_L,A,sum_A)]

    if (interventional_distribution){
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
    while (j < number_visits && nrow(people_atrisk)>0){
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
        
        ## old:         mins=Rfast::rowMins(ttt,value=FALSE)
        mins = apply(ttt,1,function(row){which.min(row)})
        people_atrisk[,event := factor(mins,levels = 1:5,labels = c("V","L","C","Y","D"))]
        ## old:         people_atrisk[,time := Rfast::rowMins(ttt,value=TRUE)]
        people_atrisk[,time := apply(ttt,1,min)]
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
        people_atrisk[event == "V",A := stats::rbinom(.N,1,propensity_A)]

        # update propensity score
        # if (interventional_distribution)
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
        ## old:     mins = Rfast::rowMins(ttt,value = FALSE)
        mins=apply(ttt,1,function(row){which.min(row)})
        p[,terminal_event := factor(mins,levels = 1:3,labels = c("C","Y","D"))]
        ## old:     p[,terminal_time := Rfast::rowMins(ttt,value = TRUE)]
        p[,terminal_time := apply(ttt,1,min)]
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
        # random baseline treatment
        # treatment_baseline <- pop[time == 0,data.table::data.table(id,date = 0, A = A_0)]
        # treatment_data <- pop[event == "V",data.table::data.table(id,date = time, A = A)]
        treatment_baseline <- pop[time == 0 & A_0 == 1,data.table::data.table(id,date = 0)]
        treatment_data <- pop[event == "V" & A == 1,data.table::data.table(id,date = time)] 
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
