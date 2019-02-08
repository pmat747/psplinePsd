#' @title Post-process of a psd object
#' @description This function allows to discard a specified number of samples as burn-in period and also thin them.
#' @param x a psd object
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (post-processing)
#' @return A list with S3 class 'psd' containing the following components:
#'    \item{psd.median,psd.mean}{psd estimates: (pointwise) posterior median and mean}
#'    \item{psd.p05,psd.p95}{90\% pointwise credibility interval}
#'    \item{psd.u05,psd.u95}{90\% uniform credibility interval}
#'    \item{fpsd.sample}{posterior power spectral density estimates}
#'    \item{k}{number of B-splines}
#'    \item{tau,phi,delta,V}{posterior traces of model parameters}
#'    \item{ll.trace}{trace of log likelihood}
#'    \item{pdgrm}{periodogram}
#'    \item{n}{integer length of input time series}
#'    \item{db.list}{posterior spectral density estimates}
#'    \item{DIC}{deviance information criterion}
#'    \item{count}{acceptance probabilities for the weigths}
#' @seealso \link{plot.psd}
#' @references Edwards, M. C., Meyer, R., and Christensen, N. (2018), Bayesian nonparametric spectral density estimation using B-spline priors, \emph{Statistics and Computing}, <https://doi.org/10.1007/s11222-017-9796-9>.
#'
#' Choudhuri, N., Ghosal, S., and Roy, A. (2004), Bayesian estimation of the spectral density of a time series, \emph{Journal of the American Statistical Association}, 99(468):1050--1059.
#'
#' @examples
#' \dontrun{
#' set.seed(123456)
#'
#' # Generate AR(1) data with rho = 0.9
#' n = 128
#' data = arima.sim(n, model = list(ar = 0.9))
#' data = data - mean(data)
#'
#' # Run MCMC (may take some time)
#' mcmc = gibbs_pspline(data, 5000, 0)
#' mcmc = burnin(mcmc, burnin = 500, thin = 10)
#' require(beyondWhittle)  # For psd_arma() function
#' freq = 2 * pi / n * (1:(n / 2 + 1) - 1)[-c(1, n / 2 + 1)]  # Remove first and last frequency
#' psd.true = psd_arma(freq, ar = 0.9, ma = numeric(0), sigma2 = 1)  # True PSD
#' plot(mcmc)  # Plot log PSD (see documentation of plot.psd)
#' lines(freq, log(psd.true), col = 2, lty = 3, lwd = 2)  # Overlay true PSD
#' }
#' @export
postprocess = function(x, burnin, thin = 1){

  N = length(x$tau);

  if(N < burnin){
    stop("burnin must be lower than the number of posterior samples");
  }

  if(class(x) != "psd"){
    stop("The object to postprocess must be a psd object");
  }

  index  = seq(from = burnin + 1, to = N, by = thin);

  cat(paste("The number of posterior samples now is ", length(index),
            sep = ""), "\n");

  fpsd.sample = x$fpsd.sample[, index];

  # Compute point estimates and 90% Pointwise CIs
  psd.median <- apply(fpsd.sample, 1, stats::median);
  psd.mean   <- apply(fpsd.sample, 1, mean);
  psd.p05    <- apply(fpsd.sample, 1, stats::quantile, probs=0.05);
  psd.p95    <- apply(fpsd.sample, 1, stats::quantile, probs=0.95);

  # Transformed versions of these for uniform CI construction
  log.fpsd.sample <- apply(fpsd.sample, 2, function(y)logfuller(y));
  log.fpsd.s      <- apply(log.fpsd.sample, 1, stats::median);
  log.fpsd.mad    <- apply(log.fpsd.sample, 1, stats::mad);
  log.fpsd.help   <- apply(log.fpsd.sample, 1, uniformmax);
  log.Cvalue      <- stats::quantile(log.fpsd.help, 0.9);

  # Compute Uniform CIs
  psd.u95 <- exp(log.fpsd.s + log.Cvalue * log.fpsd.mad);
  psd.u05 <- exp(log.fpsd.s - log.Cvalue * log.fpsd.mad);

  output = list(psd.median = psd.median,
                psd.mean   = psd.mean,
                psd.p05    = psd.p05,
                psd.p95 = psd.p95,
                psd.u05 = psd.u05,
                psd.u95 = psd.u95,
                fpsd.sample = fpsd.sample,
                k = x$k,
                tau = x$tau[index],
                phi = x$phi[index],
                delta = x$delta[index],
                V = x$V[, index],
                ll.trace = x$ll.trace[index],
                pdgrm = x$pdgrm,
                n = x$n,
                db.list = x$db.list,
                DIC = x$DIC,
                count = x$count[index]);

  class(output) = "psd"  # Assign S3 class to object

  return(output);

}

