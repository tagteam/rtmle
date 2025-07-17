### tmle_update.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:54)
## Version:
## Last-Updated: Nov 25 2024 (08:05)
##           By: Thomas Alexander Gerds
##     Update #: 47
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
tmle_update <- function(Y,
                        offset,
                        intervention_probs,
                        outcome_free_and_uncensored,
                        intervention_match,
                        Q.cont1=NULL,Q.cont0=NULL,
                        g.star=NULL,
                        intervention_type=NULL) {
  N <- length(Y)
  if (length(intervention_probs) == 0) intervention_probs <- 1
  ## FIXME: are there better ways to remove those censored in current interval?
  subjects_with_weights <- !is.na(Y) & outcome_free_and_uncensored & as.vector(intervention_match)
  weights <- numeric(N)
  if(intervention_type=="stochastic"){
    weights[subjects_with_weights] <- intervention_probs[subjects_with_weights]## the weights refers to the probability of getting intervention based on the observed data
  }
  else{
    weights[subjects_with_weights] <- 1/intervention_probs[subjects_with_weights]
  }

  if (anyNA(weights)) stop("NA in weights")
  if (any(weights > 0)) {
    f <- stats::as.formula("Y ~ -1 + S1 + offset(offset)")
    data.temp <- data.table(Y, S1 = rep(1,N),offset)
    has_weight <- weights > 0
    weights <- as.vector(scale(weights[has_weight], center = FALSE))
    m <- stats::glm(formula = f,
                    family = stats::quasibinomial(),
                    data = data.frame(data.temp[has_weight, ],weights),
                    weights = weights,
                    control = stats::glm.control(maxit = 100))

    if(intervention_type=="stochastic"){
      ## in this case, the rtmle predicted outcome is an average respect to the stochatsic intervention
      ## so we need to have the intervention at a=0 and a=1, as the prediction at a=0,a=1
      Qstar1 <- stats::predict(m, newdata = copy(data.temp)[,offset:=Q.cont1], type = "response")
      Qstar0 <- stats::predict(m, newdata = copy(data.temp)[,offset:=Q.cont0], type = "response")
      Qstar<-Qstar1*g.star + Qstar0*(1-g.star)
    }
    else{
      Qstar <- stats::predict(m, newdata = data.temp, type = "response")
    }

  } else {
    warning("No TMLE update because no subject has positive weight")
    Qstar <- stats::plogis(offset)
    if(intervention_type=="stochastic"){
    Qstar<- stats::plogis(Q.cont1)*g.star + stats::plogis(Q.cont0)*(1-g.star)
      }
    m <- "no Qstar fit because no subjects alive, uncensored, following intervention"
  }
  return(Qstar)
}



######################################################################
### tmle_update.R ends here
