### plot_adherence.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: dec 11 2025 (10:23) 
## Version: 
## Last-Updated: feb 25 2026 (16:12) 
##           By: Thomas Alexander Gerds
##     Update #: 20
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
plot_adherence <- function(x,...){
    time_nonadherence = last_interval = event_nonadherence = value = variable = first_event = event = treatment = NULL
    # find time to first deviation from regime where
    # death is a competing risk and data may be right censored
    dt_nonadherence <- do.call(rbind,lapply(names(x$protocols),function(pro){
        initiators <- x$protocols[[pro]]$intervention_match[,1,drop = TRUE] == 1
        first_deviation <- apply(x$protocols[[pro]]$intervention_match[initiators == 1,,drop = FALSE], 1, function(x) match(0, x))
        if (length(x$names$censoring)>0){
            Cvars <- grep(paste0("^",x$names$censoring,"_[0-9]+$"),names(x$prepared_data),value = TRUE)
            censored_time <- apply(x$prepared_data[initiators == 1,Cvars,with = FALSE], 1, function(x) match("censored", x))
        }
        adherence_data <- cbind(protocol = pro, data.table(x$followup[initiators == 1]),first_deviation = first_deviation,censored_time = censored_time)
        adherence_data[,time_nonadherence := pmin(last_interval,first_deviation,censored_time,na.rm = TRUE)]
        # initialize with 0
        adherence_data[,event_nonadherence := 0]
        # value 2 if not censored (competing risk or outcome)
        adherence_data[is.na(censored_time),event_nonadherence := 2]
        # value 2 when non-adherence is observed
        adherence_data[!is.na(first_deviation),event_nonadherence := 1]
        adherence_data[,data.table::data.table(protocol,time_nonadherence,event_nonadherence)]
    }))
    fit_nonadherence <- prodlim::prodlim(Hist(time_nonadherence,event_nonadherence)~protocol,
                                         data = dt_nonadherence)
    prodlim::ggprodlim(fit_nonadherence,cause = 1,ylim = c(0,100))+ ylab("Non-adherence")
}


######################################################################
### plot_adherence.R ends here
