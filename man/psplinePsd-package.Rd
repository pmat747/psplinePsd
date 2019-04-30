\name{psplinePsd-package}
\alias{psplinePsd-package}
\alias{psplinePsd}
\docType{package}
\title{
Spectral density estimation using P-spline priors
}
\description{
Implementation of a Metropolis-within-Gibbs MCMC algorithm to estimate the spectral density of a stationary time series using P-spline priors.  The algorithm makes use of the Whittle likelihood to produce pseudo-posterior samples and is based on the work presented in Maturana-Russel, P. and Meyer, R. ArXiv.  This package follows the bsplinePsd package structure.
}
\details{
The function \code{gibbs_pspline} is an implementation of the (serial version of the) MCMC algorithm presented in Maturana-Russel et al. (2019).  This algorithm uses a P-spline prior to estimate the spectral density of a stationary time series and is similar to the B-spline prior algorithm of Edwards et al. (2018).  Even though the number of knots is fixed, their location can be  efficiently determined according to the periodogram.  Alternatively, they can be equally spaced.  In addition, a pilot posterior set of samples can be used to calibrate the proposals and specify starting values.
}
\author{
Patricio Maturana-Russel, Renate Meyer

Maintainer: Patricio Maturana-Russel <p.russel@auckland.ac.nz>
}
\references{
Edwards, M. C., Meyer, R., and Christensen, N. (2018), Bayesian nonparametric spectral density estimation using B-spline priors, \emph{Statistics and Computing}, <https://doi.org/10.1007/s11222-017-9796-9>.

Maturana-Russel, P., and Meyer, R. (2019), Spectral density estimation using P-spline priors. ArXiv
}
\keyword{ package }
