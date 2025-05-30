% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{mh}
\alias{mh}
\title{This function uses the Metropolis-Hastings algorithm to sample from the posterior distributions of scaled heritability and proportion of causal genes}
\usage{
mh(
  gentres_chr_ldblock,
  lddf_chr_ldblock,
  nk,
  chain_length = 1e+05,
  burnin = 100,
  madj_likelihood = FALSE,
  tau_a = 2e-10,
  tau_x1 = 2e-09,
  tau_x2 = 2e-07,
  tau_b = 2e-04,
  beta_shape2 = 1.0001,
  verbose. = TRUE
)
}
\arguments{
\item{gentres_chr_ldblock}{gene-based association test results for a single chromosome for a random set of independent genes}

\item{lddf_chr_ldblock}{weighted LD scores for a single chromosome for a random set of independent genes}

\item{nk}{GWAS sample size}

\item{chain_length}{length of the Markov chain to use}

\item{burnin}{number of originally drawn samples to exclude from the observed posterior distribution}

\item{tau_a}{smallest parameter of the trapezoidal density function that is the prior for scaled heritability}

\item{tau_x1}{second-smallest parameter of the trapezoidal density function that is the prior for scaled heritability}

\item{tau_x2}{second-largest parameter of the trapezoidal density function that is the prior for scaled heritability}

\item{tau_b}{largest parameter of the trapezoidal density function that is the prior for scaled heritability}

\item{beta_shape2}{second shape parameter of the prior distribution of the proportion of non-causal genes}

\item{verbose.}{TRUE if results should be printed as links in the Markov chain are added, FALSE otherwise}

\item{madj}{TRUE if you want to divide likelihoods by the number of genes used, FALSE otherwise}
}
\description{
This function uses the Metropolis-Hastings algorithm to sample from the posterior distributions of scaled heritability and proportion of causal genes
}
\examples{
mh()
}
\keyword{prior}
