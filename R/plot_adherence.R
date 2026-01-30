### plot_adherence.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: dec 11 2025 (10:23) 
## Version: 
## Last-Updated: dec 12 2025 (16:07) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
plot_adherence <- function(x,...){
    # find time to first deviation from regime where
    # death is a competing risk and data may be right censored
    times <- x$time[-1]
    x$names$

    dt <- x$prepared_data[,grep(paste0("^(Lira|Placebo|Censored|primary.outcome|Dead)"),names(x$prepared_data)),with = FALSE]
    dt_time2event <- do.call(rbind,lapply(c("Lira","Placebo"),function(a){
        a_id <- which(x$prepared_data[[paste0(a,"_0")]] == 1)
        dt_a <- x$prepared_data[a_id,grep(paste0("^(",a,"|Censored|primary.outcome|Dead)"),names(x$prepared_data)),with = FALSE]
        ## need deviations to be events and hence flip the values
        A_cols <- grep(a, colnames(dt_a), value = TRUE)
        dt_a[, (A_cols) := lapply(.SD, function(x) as.numeric(x == "0")), .SDcols = A_cols]
        ## the outcome, death and censoring variables should be numeric
        C_cols <- grep("Censored", colnames(dt_a), value = TRUE)
        dt_a[, (C_cols) := lapply(.SD, function(x) as.numeric(x == "censored")), .SDcols = C_cols]
        YD_cols <- grep("primary.outcome|Dead", colnames(dt_a), value = TRUE)
        dt_a[, (YD_cols) := lapply(.SD, function(x) as.numeric(x == "1")), .SDcols = YD_cols]
        out <- dt_a[,first_event := names(.SD)[apply(.SD, 1, function(v){match(1, v)})]][,.(first_event)]
        out[,event := sapply(strsplit(first_event,"_"),"[",1)]
        out[,time := sapply(strsplit(first_event,"_"),"[",2)]
        out[,first_event := NULL]
        out[is.na(event),event := "Censored"]
        out[is.na(time),time := 11]
        out[,treatment := a]
        # change labels
        out[event == a,event := "Non-adherent"]
        out[,time:=as.numeric(as.character(factor(time,
                                                  levels=c("1","2","3","4","5","6","7","8","9","10","11"),
                                                  labels=c("1","3","6","12","18","24","30","36","42","48","54"))))]
        out
    }))
    liracolor <- "#0072B5FF"
    placebocolor <- "gray56"
    dt_time2event[,treatment:=factor(treatment,levels=c("Lira","Placebo"),labels=c("Liraglutide","Placebo"))]
    fit_nonadherence <- prodlim(Hist(time,event,cens.code = "Censored")~treatment,data = dt_time2event)
    make_panelA <- function(a){
        fit <- prodlim(Hist(time,event,cens.code = "Censored")~1,data = dt_time2event[treatment == a])
        ggprodlim(fit,cause = "Non-adherent",position_atrisk = c(0,6,12,18,24,30,36,42),ylim = c(0,40),xlim = c(0,42))+
            ylab("Non-adherence")+ xlab("Months since randomization")+
            scale_x_continuous(breaks=c(0,3,6,12,18,24,30,36,42),labels=c(0,3,6,12,18,24,30,36,42),limits=c(0,42))+
            theme_minimal(base_size = 12)+ggtitle(a)
    }
    panelA_lira <- make_panelA("Liraglutide")
    panelA_placebo <- make_panelA("Placebo")

}


######################################################################
### plot_adherence.R ends here
