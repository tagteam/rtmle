% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_effects.R
\name{compare_effects}
\alias{compare_effects}
\title{Compute estimate for varying effects}
\usage{
compare_effects(
  estimator,
  N = 1000,
  beta_L0_L = 2,
  beta_L_D = 1,
  beta_A0_L = 0,
  beta_L0_D = 0,
  nu = rep(1.1, 4),
  cens = 0,
  eta = rep(0.1, 4),
  setting = 1
)
}
\arguments{
\item{estimator}{A user specified estimator. A function that takes as input data and N,
and calculates an estimate.}

\item{N}{Number of individuals simulated}

\item{beta_L0_L}{Parameter for the effect of L0 on L.}

\item{beta_L_D}{Parameter for the effect of L on D.}

\item{beta_A0_L}{Parameter for the effect of L0 on L.}

\item{beta_L0_D}{Parameter for the effect of L0 on D.}

\item{nu}{Vector of length 4 of scale parameters for the Weibull intensity with parameterization
\deqn{\eta \nu t^{\nu - 1}}. Default is set to 1.1 for all events.}

\item{cens}{Indicator for whether censoring is present}

\item{eta}{Vector of length 4 of shape parameters for the Weibull intensity. Default is set to
0.1 for all events.}

\item{setting}{Setting specifier}
}
\value{
Plot
}
\description{
`compare_effects` is a function for computing a user specified estimate for three
different effects of A0 on the intensity of L, and for a varying parameter.
}
\examples{
# First we define an estimator
estimator1 <- function(data, N) {

T2D_events <- data[Delta == 3]
return(c(nrow(T2D_events[A0 == 1 & Time < 5])/N, nrow(T2D_events[A0 == 0 & Time < 5])/N))

}

res1 <- compare_effects(estimator = estimator1, N = 100, beta_L_D = seq(0,1,by = 0.1))
}
