#' @title Post-process of a psd object
#' @description This function allows to discard a specified number of samples as burn-in period and also thin them.
#' @param x a psd object
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (post-processing)
#' @return A list with S3 class 'psd' containing the following updated components:
#'    \item{psd.median,psd.mean}{psd estimates: (pointwise) posterior median and mean}
#'    \item{psd.p05,psd.p95}{90\% pointwise credibility interval}
#'    \item{psd.u05,psd.u95}{90\% uniform credibility interval}
#'    \item{fpsd.sample}{posterior power spectral density estimates}
#'    \item{tau,phi,delta,V}{posterior traces of model parameters}
#'    \item{ll.trace}{trace of log likelihood}
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
#' n = 128;
#' data = arima.sim(n, model = list(ar = 0.9));
#' data = data - mean(data);
#'
#' # Run MCMC (may take some time)
#' mcmc = gibbs_pspline(data, 5000, 0);
#' mcmc = burnin(mcmc, burnin = 500, thin = 10);
#'
#' require(beyondWhittle)  # For psd_arma() function
#' freq = 2 * pi / n * (1:(n / 2 + 1) - 1)[-c(1, n / 2 + 1)]  # Remove first and last frequency
#' psd.true = psd_arma(freq, ar = 0.9, ma = numeric(0), sigma2 = 1)  # True PSD
#' plot(mcmc)  # Plot log PSD (see documentation of plot.psd)
#' lines(freq, log(psd.true), col = 2, lty = 3, lwd = 2)  # Overlay true PSD
#'
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

  # parameters

  tau = x$tau[index];
  V   = x$V[, index];

  ###########
  ### DIC ###
  ###########

  # anSpecif$FZ: fast_ft(data);# FFT data to frequency domain. NOTE: Must be mean-centred.
  FZ    <- x$anSpecif$FZ;
  pdgrm <- abs(FZ) ^ 2; # Periodogram: NOTE: the length is n here.
  omega <- 2 * (1:(x$anSpecif$n / 2 + 1) - 1) / x$anSpecif$n;  # Frequencies on unit interva

  tau_mean = mean(tau);

  v_means = unname(apply(V, 1, mean));

  l = llike(omega, FZ, x$anSpecif$k, v = v_means,
            tau = tau_mean, pdgrm, x$anSpecif$degree, x$db.list)$llike;

  ls = apply(rbind(tau, V), 2, function(y){
             llike(omega, FZ, x$anSpecif$k, v = y[-1],
             tau = y[1], pdgrm, x$anSpecif$degree, x$db.list)$llike});
  ls = unname(ls);

  # http://kylehardman.com/BlogPosts/View/6
  # DIC = -2 * (l - (2 * (l - mean(ls))));

  D_PostMean = -2 * l;
  D_bar      = -2 * mean(ls);
  pd         = D_bar - D_PostMean;

  DIC = list(DIC = 2 * D_bar - D_PostMean, pd = pd);

  x[["psd.median"]] = psd.median;
  x[["psd.mean"]]   = psd.mean;
  x[["psd.p05"]]    = psd.p05;
  x[["psd.p95"]]    = psd.p95;
  x[["psd.u05"]]    = psd.u05;
  x[["psd.u95"]]    = psd.u95;
  x[["fpsd.sample"]]= fpsd.sample;
  #x[["k"]]          = x$k;
  x[["tau"]]        = tau; # defined above
  x[["phi"]]        = x$phi[index];
  x[["delta"]]      = x$delta[index];
  x[["V"]]          = V; # defined above
  x[["ll.trace"]]   = x$ll.trace[index];
  #x[["pdgrm"]]      = x$pdgrm;
  #x[["n"]]          = x$n;
  #x[["db.list"]]    = x$db.list;
  x[["DIC"]]        = DIC;
  x[["count"]]      = x$count[index];

  return(x);

}

