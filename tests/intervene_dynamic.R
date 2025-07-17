## example for dynamic rule for the intervene_function,
## I believe that it is only needed to have the value at the actual time
## because the weight will be non-zero only if intervention.match=1
## meaning for those individuals that actually follow the regimen


### this is going to change only the value for the current time
### not the one before or after


## we could make define this externally, giving an example
## and this can depends either on the time we are at the moment
## or on all the previous ones

intervention_fun<-function(db,treatment,time){
  ## we can specify, if treatment is A--->
  ## if it is B---> other definition
  # if(time==0){
  #   set.seed(389)
  #   values<-rbinom(nrow(db),1,0.5)
  # }
  # else{
  # #cov<-paste0("age_",time)
  #   cov<-"age"
  # values<-ifelse(db[,cov,with = FALSE]> 70,1,0)}
  values=rep(1, nrow(db))
  factor(values, levels = c("0","1"))
  }


### like  this we intervene on all treatment at once, for each time!
intervene_dynamic <- function(data,
                      intervention_table){

  interdata <- copy(data)
  N <- NROW(interdata)
  for (k in 1:nrow(intervention_table)){
    set(interdata,
        j = intervention_table[k][["variable"]],
        value = intervention_fun(db=interdata,
                                 treatment=intervention_table[k][["variable"]],
                                 time=intervention_table[k][["time"]]))
  }


  interdata

}
