### intervene_stochastic_SI.R ---
#----------------------------------------------------------------------
## Author: Alessandra Meddis
## Created: Nov 27 2024 (15:30)
## Version:
## Last-Updated: 11 Dec 2024 (08:59)
##           By: Alessandra
##     Update #: 155
#----------------------------------------------------------------------
##
### Commentary: when assigning a stochastic intervention, we need the definition of
#               the probability distribution for the Treatment assignment
##
### Change Log:
#----------------------------------------------------------------------


## we provide directly the trt.name because in case of multiple treatment, only the current.time is not useful
intervene_stochastic<-function(data, name.treat, current.time){


  ## in case of more then one treatment, it is needed!
  # prob is the probability pf receiving treatment if this thing happenes
  ## we can also have prob<-0.5 so it is random every time regardless of everything!
  # if(grep("A",name.trt)){
  #
  #   age_var<-paste0("age_",current.time,"", sep="")
  #   past_trt<-paste0("A_",current.time-1,"", sep="")
  #   # if they are under treatment
  #   prob<-ifelse(data[,..past_trt]==0, ifelse(data[, ..age_var]>30, 0.2,0.8), 1)
  #
  # }

  ### define the treatment based on the age and on the past treatment
  ### we would like to say that if a person had received treatment then she stays on
  past_trt<-name.treat
  if(current.time>0){
    past_trt<-paste0("A_",current.time-1,"", sep="")
  }
  prob<-ifelse(data[,..past_trt]==0, ifelse(data$age>30, 0.4,0.6), 1)
  #prob=rep(1,nrow(data))
  prob
}
