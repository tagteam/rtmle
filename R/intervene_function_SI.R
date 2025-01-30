### intervene_function_SI.R ---
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
#               this is an example for the definition of the stochastic intervention that
#               needs to be specified in the protocol function
##
### Change Log:
#----------------------------------------------------------------------


# in data we have the current_data--> the data that we have at this specific time
# current.time is the time considered at this point

intervene_function<-function(data,intervention_table,current.time){
  if(current.time==1)
  { value<-rep(1,nrow(data))
  }
  else
  {
    name.treat<-intervention_table[time==current.time-1, ]$variable
    value<-ifelse(data[,..name.treat]==0,0, 0.8) # once they stop treatment, they cannot re-start it otherwise they have 80% prb to stay under tretament
  }
  #name.treat<-intervention_table[time==current.time, ]$variable
  # we would like to set this to 1 if under the specified treatment and to 0 otehrwise (as for the static intervention)
  # check the variable:
  #value<-ifelse(data[,..name.treat]==1,1, 0)
  #value
  rep(0,nrow(data))
}

