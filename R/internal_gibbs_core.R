#' Generate a B-spline density basis of any degree
#'
#' "dbspline" function has been taken from "bsplinePsd" package
#'    Other functions have been modified for our case.
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

  expV   <- exp(v)
  #weight <- (diag(k-1) - expW %*% t(rep(1, k-1))/(1+sum(expV)) ) %*% expV;
  weight <- expV / (1+sum(expV));
  weight <- c(weight, 1-sum(weight));

  psd <- densityMixture(weight, db.list)
  epsilon <- 1e-20
  psd <- pmax(psd, epsilon)
  psd <- psd[-c(1, length(psd))]
  return(list(psd = psd, db.list = db.list))
}

#' This function can deal with Inf values when taking exponential of weigths,
#'   unlike the function "qpsd" defined above
#' @importFrom Rcpp evalCpp
#' @keywords internal
qpsd = function (omega, k, v, degree, db.list)
{
  v      <- as.numeric(v);
  expV   <- exp(v);

  # converting v into a weight

  if(any(is.infinite(expV))){

    # using a log-version in case w contains large values

    ls     = logplusvec(c(0,v)); # log(1+ sum e^{v_i})
    weight = exp(v - ls);

  }else{
    # Matrix form
    # weight <- (diag(k-1) - expV %*% t(rep(1, k-1))/(1+sum(expV)) ) %*% expV;
    # weight <- c(weight, 1-sum(weight));
    ls     = 1 + sum(expV);
    weight = expV / ls;
  }
  s       <- 1-sum(weight);
  weight  <- c(weight, ifelse(s<0, 0, s));

  psd     <- densityMixture(weight, db.list);
  epsilon <- 1e-20;
  psd     <- pmax(psd, epsilon);
  psd     <- psd[-c(1, length(psd))];
  return(list(psd = psd, db.list = db.list));
}

#' @keywords internal
logplus <- function(x,y)
{
  if(x>y) x + log(1+exp(y-x))
  else    y + log(1+exp(x-y))
}

# vec-version
#' @keywords internal
logplusvec = function(x){
  r = -Inf;
  for(i in x){
    r = logplus(r, i);
  }
  return(r);
}

#' log-prior
#' @keywords internal
lprior = function (k, v, tau, tau.alpha, tau.beta, phi, phi.alpha, phi.beta,
                   delta, delta.alpha, delta.beta, P)
{
  # Sigma^(-1) = P

  logprior <- (k-1)*log(phi)/2 - phi * t(v) %*% P %*% v / 2 + #MNormal on weights

    phi.alpha * log(delta) + (phi.alpha - 1) * log(phi) - phi.beta * delta * phi + # log prior for phi (Gamma distb)

    (delta.alpha - 1) * log(delta) - delta.beta * delta - # log prior for delta

    (tau.alpha + 1) * log(tau) - tau.beta/tau; # log prior for tau (Inverse Gamma)

  return(logprior)
}

#' log Whittle likelihood
#' @keywords internal
llike = function (omega, FZ, k, v, tau, pdgrm, degree, db.list)
{
  n <- length(FZ)
  m <- n - 2
  qq.psd <- qpsd(omega, k, v, degree, db.list)
  q.psd <- qq.psd$psd
  q <- rep(NA, m)
  q[1] <- q.psd[1]
  q[m] <- q.psd[length(q.psd)]
  q[2 * 1:(m/2 - 1)] <- q[2 * 1:(m/2 - 1) + 1] <- q.psd[1:(m/2 -
                                                             1) + 1]
  f <- tau * q
  llike <- -sum(log(f) + pdgrm[2:(n - 1)]/(f * 2 * pi))/2
  return(llike)
}

#' Unnormalised log posterior
#' @keywords internal
lpost = function (omega, FZ, k, v, tau, tau.alpha, tau.beta,
                  phi, phi.alpha, phi.beta, delta, delta.alpha, delta.beta,
                  P, pdgrm, degree, db.list)
{
  ll <- llike(omega, FZ, k, v, tau, pdgrm, degree, db.list)

  lp <- ll + lprior(k, v, tau, tau.alpha, tau.beta, phi,
                    phi.alpha, phi.beta, delta, delta.alpha, delta.beta, P)

  return(lp)
}
