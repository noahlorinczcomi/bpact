% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{ptrap}
\alias{ptrap}
\title{This function calculates the cumulative density of the trapezoidal distribution}
\usage{
ptrap(x, a, x1, x2, b, log = FALSE)
}
\arguments{
\item{x}{argument to the cumulative trapezoidal density function}

\item{a}{smallest parameter of the trapezoidal density function}

\item{x1}{second-smallest parameter of the trapezoidal density function}

\item{x2}{second-largest parameter of the trapezoidal density function}

\item{b}{largest parameter of the trapezoidal density function}

\item{log}{TRUE if log density should be returned, FALSE otherwise}
}
\description{
This function calculates the cumulative density of the trapezoidal distribution
}
\examples{
ptrap(2,0,1,2,3)
ptrap()
}
\keyword{trapezoidal}
