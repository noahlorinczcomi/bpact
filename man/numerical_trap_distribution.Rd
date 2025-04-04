% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{numerical_trap_distribution}
\alias{numerical_trap_distribution}
\title{This function numerically estimates the parameters of a trapezoidal distribution using observations}
\usage{
numerical_trap_distribution(
  taus,
  subn = length(taus),
  atrim = 0,
  btrim = 0.975,
  doplot = FALSE
)
}
\arguments{
\item{taus}{vector of observed random variables following a trapezoidal distribution}

\item{subn}{number of sub-sampled observations to use if length of taus is very large}

\item{atrim}{the lower quantile of of observed values to consider, such as 0.05}

\item{btrim}{the upper quantile of of observed values to consider, such as 0.95}

\item{doplot}{TRUE if a plot should be generated, FALSE otherwise}
}
\description{
This function numerically estimates the parameters of a trapezoidal distribution using observations
}
\examples{
numerical_trap_distribution(c(1,3,4,5,7))
numerical_trap_distribution()
}
\keyword{trapezoidal}
