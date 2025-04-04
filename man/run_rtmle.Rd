% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_rtmle.R
\name{run_rtmle}
\alias{run_rtmle}
\title{Sequential regression with TMLE update step for discretized follow-up data}
\usage{
run_rtmle(
  x,
  targets,
  learner = "learn_glm",
  time_horizon,
  refit = FALSE,
  seed = NULL,
  ...
)
}
\arguments{
\item{x}{object of class \code{rtmle}}

\item{targets}{Selection of targets to be analysed. If missing all targets in x$targets are analysed.}

\item{learner}{A function which is called to fit (learn) the nuisance parameter models.}

\item{time_horizon}{The time horizon at which to calculate risks. If it is a vector the analysis will be performed for each element of the vector.}

\item{refit}{Logical. If \code{TRUE} ignore any propensity score and censoring models learned in previous calls to this function. Default is \code{FALSE}.}

\item{seed}{Seed used for cross-fitting}

\item{...}{Additional arguments passed to the learner function.}
}
\value{
The modified object contains the fitted nuisance parameter models and the estimate of the target parameter.
}
\description{
This function runs the analysis defined in previous steps.
}
\examples{
# ------------------------------------------------------------------------------------------
# Intervening on a single treatment variable
# ------------------------------------------------------------------------------------------

set.seed(17)
tau <- 3
ld <- simulate_long_data(n = 91,number_visits = 20,
                         beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),
                         register_format = TRUE)
x <- rtmle_init(intervals = tau,name_id = "id",name_outcome = "Y",name_competing = "Dead",
                name_censoring = "Censored",censored_label = "censored")
x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
protocol(x) <- list(name = "Always_A",
                    intervention = data.frame("A" = factor("1",levels = c("0","1"))))
protocol(x) <- list(name = "Never_A",
                    intervention = data.frame("A" = factor("0",levels = c("0","1"))))
prepare_data(x) <- list()
target(x) <- list(name = "Outcome_risk",
                  strategy = "additive",
                  estimator = "tmle",
                  protocols = c("Always_A","Never_A"))
x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau)
summary(x)



}
\author{
Thomas A Gerds \email{tag@biostat.ku.dk}
}
