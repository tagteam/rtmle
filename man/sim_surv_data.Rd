% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_surv_data.R
\name{sim_surv_data}
\alias{sim_surv_data}
\title{Simulate Survival Data}
\usage{
sim_surv_data(N, beta = NULL, eta = rep(0.1, 2), nu = rep(1.1, 2))
}
\arguments{
\item{N}{A double of the number of individuals}

\item{beta}{A 2X2 matrix with the effect of L0 and A0 on Censoring and Event.
The columns represent Censoring and Death, while the rows represent L0 and A0.}

\item{eta}{Vector of  length 2 of shape parameters for the Weibull hazard with parameterization
\deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.}

\item{nu}{Vector of length 2 of scale parameters for the Weibull hazard.}
}
\value{
Data frame containing the simulated survival data
}
\description{
Function to simulate data from a Survival setting. The function simulates data corresponding to $N$ individuals that
are at risk for censoring (0) and an event (1).
}
\examples{
sim_surv_data(10)
}
