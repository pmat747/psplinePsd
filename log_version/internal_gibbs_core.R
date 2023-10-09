#' Generate a B-spline density basis of any degree
#' "dbspline" function has been taken from "bsplinePsd" package
#' The other functions have been modified for our case.
#'
#' @importFrom splines splineDesign
#' @importFrom Rcpp evalCpp
#' @keywords internal
dbspline = function (x, knots, degree = 3)  {

  knots.mult <- c(rep(knots[1], degree), knots, rep(knots[length(knots)], degree))
  nknots = length(knots.mult)  # Number of knots including external knots

  B <- splines::splineDesign(knots.mult, x, ord = degree + 1, outer.ok = TRUE)

  # Trivial normalisation formula
  bs_int = (knots.mult[-(1:(degree + 1))] - knots.mult[-((nknots - degree):nknots)]) / (degree + 1)
  if (any(bs_int == 0)) bs_int[which(bs_int == 0)] = Inf  # Makes B.norm = 0 rather than NaN

  B.norm <- t(B) / bs_int  # Normalise

  return(B.norm)

}

#' Compute unnormalised PSD using random mixture of B-splines
#' @importFrom Rcpp evalCpp
#' @keywords internal
qpsd = function (omega, k, v, degree, db.list)
{
  
  psd <- densityMixture(v, db.list)
  
  return(psd)
}

#' log-prior
#' @keywords internal
lprior = function (k, v, tau, tau.alpha, tau.beta, phi, phi.alpha, phi.beta,
                   delta, delta.alpha, delta.beta, P)
{
  # Sigma^(-1) = P

  logprior <- k * log(phi)/2 - phi * t(v) %*% P %*% v / 2 + #MNormal on weights

    dgamma(phi, phi.alpha, delta * phi.beta, log = TRUE) + # log prior for phi
    
    dgamma(delta, delta.alpha, delta.beta, log = TRUE) + # log prior for delta

    dnorm(tau, 0, 100, log = TRUE); # prior for tau

  return(logprior)
}

#' log Whittle likelihood
#' @importFrom Rcpp evalCpp
#' @useDynLib psplinePsd, .registration = TRUE
#' @keywords internal
llike = function (omega, FZ, k, v, tau, tau.alpha, pdgrm, degree, db.list,
                  spec_ar)
{
  n <- length(FZ);

  # Which boundary frequencies to remove from likelihood computation
  if (n %% 2) {  # Odd length time series
    bFreq <- 1  # Remove first
  }
  else {  # Even length time series
    bFreq <- c(1, n)  # Remove first and last
  }

  # Un-normalised PSD (defined on [0, 1])
  qq.psd <- qpsd(omega, k, v, degree, db.list);

  q = unrollPsd(qq.psd, n); # Unrolls the unnormalised PSD to length n

  # Normalised PSD (defined on [0, pi])
  f <- q + tau; # CHANGED
  f <- f + spec_ar; # Correction
  
  # Whittle log-likelihood
  #llike <- -sum(log(f[-bFreq]) + pdgrm[-bFreq] / (f[-bFreq] * 2 * pi)) / 2; # original
  #llike <- -sum(f[-bFreq] + pdgrm[-bFreq] / (exp(f[-bFreq]) * 2 * pi)) / 2;
  llike <- -sum(f[-bFreq] + 
                exp(log(pdgrm[-bFreq]/(2 * pi)) - f[-bFreq]) ) / 2;
  
  return(llike)

}

#' Unnormalised log posterior
#' @keywords internal
lpost = function (omega, FZ, k, v, tau, tau.alpha, tau.beta,
                  phi, phi.alpha, phi.beta, delta, delta.alpha, delta.beta,
                  P, pdgrm, degree, db.list,
                  spec_ar)
{
  
  ll  <- llike(omega, FZ, k, v, tau, tau.alpha, pdgrm, degree, db.list, spec_ar)

  lpr <- lprior(k, v, tau, tau.alpha, tau.beta, phi,
                phi.alpha, phi.beta, delta, delta.alpha, delta.beta, P)

  #print(ll)
  #print(lpr)

  lp <- ll + lpr;

  #return(lp)
  return(list(ll = ll, lp = lp));
}
