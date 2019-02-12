#' @title Metropolis-within-Gibbs sampler for spectral inference of a stationary time series using p-splines
#' @description This function uses the Whittle likelihood and obtains samples from the pseudo-posterior to infer the spectral density of a stationary time series.
#' @details The function \code{gibbs_pspline} is an implementation of the (serial version of the) MCMC algorithm presented in Edwards et al. (2018).  This algorithm uses P-splines to estimate the spectral density of a stationary time series and can be considered a particular case of the algorithm of Edwards et al. (2018), which used a B-spline prior allowing the number of B-splines functions to be variable.
#'          We consider the prior on the spectral density given by
#'          \deqn{f(w) = \tau \sum_{j=1}^{k}w_{j}B_{j}(w)}
#'          where \eqn{B_{j}} is the B-spline density.  The following prior is allocated indirectly on the weights \eqn{w}:
#'
#'          \deqn{v|\phi \delta ~ N_{k-1}(0, (\phi D^\top D)^{-1})}
#'
#'          \deqn{\phi|\delta ~ Gamma(\alpha_{\phi}, \delta \beta_{\phi})}
#'
#'          \deqn{\delta ~ Gamma(\alpha_{\delta}, \beta_{\delta})}
#'
#'          where \deqn{v_{j} = \log ( \frac{w_{j}}{1-\sum_{j=1}^{k-1} w_{j}} )}.
#' @param data numeric vector
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (post-processing)
#' @param tau.alpha,tau.beta prior parameters for tau (Inverse-Gamma)
#' @param phi.alpha,phi.beta prior parameters for phi (Gamma)
#' @param delta.alpha,delta.beta prior parameters for delta (Gamma), which is a factor of the shape hyperparameter in the Gamma prior for tau
#' @param k number of B-splines
#' @param degree positive integer specifying the degree of the B-spline densities (default is 3)
#' @param diffMatrixOrder  order of the difference penalty matrix (1, 2 or 3)
#' @param psd output from \code{gibbs_pspline} function
#' @param add logical value indicating wether to add pilot posterior samples "psd" to the current analysis
#' @return A list with S3 class 'psd' containing the following components:
#'    \item{psd.median,psd.mean}{psd estimates: (pointwise) posterior median and mean}
#'    \item{psd.p05,psd.p95}{90\% pointwise credibility interval}
#'    \item{psd.u05,psd.u95}{90\% uniform credibility interval}
#'    \item{fpsd.sample}{posterior spectral density estimates}
#'    \item{anSpecif}{a list with some of the specifications of the analysis}
#'    \item{n}{integer length of input time series}
#'    \item{tau,phi,delta,V}{posterior traces of model parameters}
#'    \item{ll.trace}{trace of log likelihood}
#'    \item{pdgrm}{periodogram}
#'    \item{db.list}{B-spline densities}
#'    \item{DIC}{deviance information criterion}
#'    \item{count}{acceptance probabilities for the weigths}
#' @seealso \link{plot.psd}
#' @references Edwards, M. C., Meyer, R., and Christensen, N. (2018), Bayesian nonparametric spectral density estimation using B-spline priors, \emph{Statistics and Computing}, <https://doi.org/10.1007/s11222-017-9796-9>.
#'
#' Choudhuri, N., Ghosal, S., and Roy, A. (2004), Bayesian estimation of the spectral density of a time series, \emph{Journal of the American Statistical Association}, 99(468):1050--1059.
#'
#' @examples
#' \dontrun{
#'
#' set.seed(1)
#'
#' # Generate AR(1) data with rho = 0.9
#' n = 128
#' data = arima.sim(n, model = list(ar = 0.9))
#' data = data - mean(data)
#'
#' # Run MCMC (may take some time)#'
#' pilotmcmc = gibbs_pspline(data, 2500, 500); # pilot run
#' mcmc1 = gibbs_pspline(data, 3000, 2000, psd = pilotmcmc);
#' mcmc2 = gibbs_pspline(data, 3000, 0, psd = mcmc1, add = TRUE); # reciclying mcmc1 samples
#'
#' require(beyondWhittle)  # For psd_arma() function
#' freq = 2 * pi / n * (1:(n / 2 + 1) - 1)[-c(1, n / 2 + 1)]  # Remove first and last frequency
#' psd.true = psd_arma(freq, ar = 0.9, ma = numeric(0), sigma2 = 1)  # True PSD
#'
#' plot(mcmc1)  # Plot log PSD (see documentation of plot.psd)
#' lines(freq, log(psd.true), col = 2, lty = 3, lwd = 2)  # Overlay true PSD
#'
#' plot(mcmc2)  # Plot log PSD (see documentation of plot.psd)
#' lines(freq, log(psd.true), col = 2, lty = 3, lwd = 2)  # Overlay true PSD
#'
#' }
#' @importFrom Rcpp evalCpp
#' @useDynLib psplinePsd, .registration = TRUE
#' @export
gibbs_pspline <- function(data,
                          Ntotal,
                          burnin,
                          thin = 1,
                          tau.alpha = 0.001,
                          tau.beta = 0.001,
                          phi.alpha = 1,
                          phi.beta = 1e-5,
                          delta.alpha = 1,
                          delta.beta = 1,
                          k = NULL,
                          degree = 3,
                          diffMatrixOrder = 3,
                          psd = NULL,
                          add = FALSE) {

  if(is.null(psd)){

    out = gibbs_pspline_simple(data = data,
                               Ntotal = Ntotal,
                               burnin = burnin,
                               thin = thin,
                               tau.alpha = tau.alpha,
                               tau.beta = tau.beta,
                               phi.alpha = phi.alpha,
                               phi.beta = phi.beta,
                               delta.alpha = delta.alpha,
                               delta.beta = delta.beta,
                               k = k,
                               degree = degree,
                               diffMatrixOrder = diffMatrixOrder);

  }else{

    out = gibbs_pspline_postProposal(data = data,
                                     Ntotal = Ntotal,
                                     burnin= burnin,
                                     thin = thin,
                                     tau.alpha = tau.alpha,
                                     tau.beta = tau.beta,
                                     phi.alpha = phi.alpha,
                                     phi.beta = phi.beta,
                                     delta.alpha = delta.alpha,
                                     delta.beta = delta.beta,
                                     k = k,
                                     degree = degree,
                                     diffMatrixOrder = diffMatrixOrder,
                                     psd = psd,
                                     add = add);
  }

  return(out);

}  # Close function
