% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gibbs_pspline_postProposal.R
\name{gibbs_pspline_postProposal}
\alias{gibbs_pspline_postProposal}
\title{This function takes into account a posterior pilot sample to calibrate
  the proposals for the weights}
\usage{
gibbs_pspline_postProposal(
  data,
  Ntotal,
  burnin,
  thin = 1,
  tau.alpha,
  tau.beta,
  phi.alpha,
  phi.beta,
  delta.alpha,
  delta.beta,
  k,
  eqSpacedKnots,
  degree,
  diffMatrixOrder,
  printIter = 100,
  psd,
  add = FALSE
)
}
\description{
This function takes into account a posterior pilot sample to calibrate
  the proposals for the weights
}
\keyword{internal}
