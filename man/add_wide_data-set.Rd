% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_wide_data.R
\name{add_wide_data<-}
\alias{add_wide_data<-}
\title{Adding wide format data}
\usage{
add_wide_data(x, ...) <- value
}
\arguments{
\item{x}{object of class \code{rtmle}}

\item{...}{Not used (not yet)}

\item{value}{Named list of data.frames (or data.tables or tibbles). Two possible names 'outcome_data' and 'timevar_data'.}
}
\value{
The modified object.
}
\description{
Adding wide format data to a rtmle object
}
\details{
This function adds a list of datasets in wide format (one line per subject) to an existing rtmle object
}
\examples{
## FIXME
}
\seealso{
\link[rtmle]{add_baseline_data<-}, \link[rtmle]{add_long_data<-}
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
