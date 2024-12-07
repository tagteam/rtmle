% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_comp_risk_data.R
\name{sim_comp_risk_data}
\alias{sim_comp_risk_data}
\title{Simulate Competing risk Data}
\usage{
sim_comp_risk_data(N, beta = NULL, eta = rep(0.1, 3), nu = rep(1.1, 3))
}
\arguments{
\item{N}{A double of the number of individuals}

\item{beta}{A 2X3 matrix with the effect of L0 and Treatment on two competing processes and censoring.
The columns represent Process 1, Process 2 and Censoring, while the rows represent L0 and Treatment.}

\item{eta}{Vector of  length 3 of shape parameters for the Weibull hazard with parameterization
\deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.}

\item{nu}{Vector of length 3 of scale parameters for the Weibull hazard.}
}
\value{
Data frame containing the simulated competing risk data
}
\description{
Function to simulate Competing Risk data.
}
\examples{
sim_comp_risk_data(10)
}
