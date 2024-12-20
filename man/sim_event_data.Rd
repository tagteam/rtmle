% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_event_data.R
\name{sim_event_data}
\alias{sim_event_data}
\title{Simulate Event Data}
\usage{
sim_event_data(
  N,
  beta = NULL,
  eta = rep(0.1, 4),
  nu = rep(1.1, 4),
  at_risk = NULL,
  term_deltas = c(1, 2)
)
}
\arguments{
\item{N}{A double for the number of simulated individuals}

\item{beta}{A matrix of doubles for the effects on the intensities. The columns represent the events. In the
default case Operation, Death, Censoring, and Covariate Change. The rows represent the baseline covariate \eqn{L0},
the indicator for operation \eqn{A}, baseline treatment \eqn{A0}, and the indicator for change in the covariate
process \eqn{L}.The \eqn{beta} matrix is by default set to 0.}

\item{eta}{Vector of shape parameters for the Weibull intensity with parameterization
\deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.}

\item{nu}{Vector of scale parameters for the Weibull intensity. Default is set to 0.1 for all events.}

\item{at_risk}{At risk function. Default is set to the operation setting. A survival or competing risk
setting can be specified as well. The \code{at_risk} function is an indicator of whether an individual is at risk
for a specific event. The function takes as input i (the index belonging to a particular individual) L (indicator
for a change in covariate) and A (indicator for change in treatment process).The function returns a vector of 0's
and 1's corresponding to whether individual i is at risk for a particular event.}

\item{term_deltas}{Terminal events. Default is set so that event 1 and 2 are terminal events.}
}
\value{
Data frame containing the simulated data. There is a column for ID, time of event (Time),
event type (Delta), baseline covariate (L0), additional covariate (L), Baseline Treatment (A0)
and Treatment Process (A).
}
\description{
\code{sim_event_data} is a function to simulate event data, e.g. observational healthcare data. The number of events
simulated corresponds to the length of the \eqn{eta} and \eqn{nu} vector, and the number of columns in the
\eqn{beta} matrix. By default 4 different types of events are simulated, chosen to represent Operation (0),
Death (1), Censoring(2) and Change in Covariate Process(3). Death and Censoring are terminal events and Operation
and Change in Covariate Process can occur once. The intensities of the various events depend upon previous events
and the pre specified \eqn{beta}, \eqn{nu} and \eqn{eta} parameters.
}
\details{
Different settings can be specified: Survival setting, competing risk setting and operation setting.
The different settings can be specified by the at_risk function. Default is the operation setting.
}
\examples{
sim_event_data(N = 10)
}
