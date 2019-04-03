#' This function does not need a posterior pilot sample
#'
#' @importFrom Rcpp evalCpp
#' @importFrom stats median rnorm runif rgamma
#' @useDynLib psplinePsd, .registration = TRUE
#' @keywords internal
gibbs_pspline_simple <- function(data,
                                 Ntotal,
                                 burnin,
                                 thin = 1,
                                 tau.alpha = 0.001,
                                 tau.beta = 0.001,
                                 phi.alpha = 1,
                                 phi.beta = 1,
                                 delta.alpha = 1e-04,
                                 delta.beta = 1e-04,
                                 k = NULL,
                                 eqSpacedKnots = FALSE,
                                 degree = 3,
                                 diffMatrixOrder = 2,
                                 printIter = 100) {

  if (burnin >= Ntotal) stop("burnin must be less than Ntotal")
  if (any(c(Ntotal, burnin, thin) %% 1 != 0 || any(c(Ntotal, burnin, thin) < 0)))
    if (any(c(tau.alpha, tau.beta) <= 0)) stop(" tau.alpha and tau.beta must be strictly positive")
  if (any(c(Ntotal, thin) %% 1 != 0) || any(c(Ntotal, thin) <= 0)) stop("Ntotal must be strictly positive integers")
  if ((burnin %% 1 != 0) || (burnin < 0)) stop("burnin must be a non-negative integer")

  n <- length(data);
  if (n %% 2 != 0) stop("this version of bsplinePsd must have n even")

  if(is.null(k)){
    k = min(round(n/4), 40);
    cat(paste("Number of B-splines k=", k, sep=""), "\n");
  }

  if( (Ntotal - burnin)/thin < k){stop("Change corresponding specifications in order to have (Ntotal-burnin)/thin > k")}

  if( (printIter<=0) || (printIter %% 1 != 0) )stop("printIter must be a positive integer value");

  # Tolerance for mean centering
  tol <- 1e-4;

  # Mean center
  if (abs(mean(data)) > tol) {
    data <- data - mean(data);
    warning("Data has been mean-centred")
  }

  # Optimal rescaling to prevent numerical issues
  rescale = stats::sd(data);
  data = data / rescale;  # Data now has standard deviation 1

  FZ <- fast_ft(data); # FFT data to frequency domain.  NOTE: Must be mean-centred.

  pdgrm <- abs(FZ) ^ 2; # Periodogram: NOTE: the length is n here.

  omega <- 2 * (1:(n / 2 + 1) - 1) / n;  # Frequencies on unit interval
  lambda <- pi * omega; # Angular frequencies on [0, pi]

  # number of samples to be stored
  N = round(Ntotal/thin);

  # Open objects for storage
  tau   <- rep(NA, N);
  phi   <- rep(NA, N); # scalar in Multivariate Normal
  delta <- rep(NA, N);

  # Difference Matrix
  P       = diffMatrix(k-1, d = diffMatrixOrder); # Third order penalty
  P       = t(P) %*% P;
  epsilon = 1e-6; #1e-6;
  P       = P + epsilon * diag(dim(P)[2]); # P^(-1)=Sigma (Covariance matrix)

  # Starting values
  tau[1]   <- stats::var(data) / (2 * pi);
  delta[1] <- delta.alpha/delta.beta;
  phi[1]   <- phi.alpha/(phi.beta * delta[1]);

  # starting value for the weights
  #w      <- rep(0,k-1);
  w = pdgrm / sum(pdgrm);
  #w = w[round(seq(1, length(w), length = k-1))];
  w = w[round(seq(1, length(w), length = k))];
  w[which(w==0)] = 1e-50; # prevents errors when there are zeros
  w = w/sum(w);
  w = w[-k];
  v = log(w / (1 - sum(w)));
  V = matrix(v, ncol = 1);

  ###
  # Since the number of knots are fixed,
  #  the B-splines only need to be calculated once.

  #newk    <- k - degree + 1;
  #knots   <- seq(0,1, length = newk);
  knots = knotLoc(data = data, k = k, degree = degree, eqSpaced = eqSpacedKnots);

  db.list <- dbspline(omega, knots, degree);
  ###

  # Store log likelihood
  ll.trace <- NULL;

  Count = NULL; # ACCEPTANCE PROBABILITY
  sigma = 1;    # variance of proposal distb for weights
  count = 0.4;  # starting value for acc pbb - optimal value
  k1    = k - 1;

  # Random values
  Zs = stats::rnorm((Ntotal-1)*k1);
  Zs = matrix(Zs, nrow = Ntotal-1, ncol = k1);
  Us = log(stats::runif((Ntotal-1)*k1, 0, 1));
  Us = matrix(Us, nrow = Ntotal-1, ncol = k1);

  # Initial values for the proposals
  phi.store   = phi[1];
  tau.store   = tau[1];
  delta.store = delta[1];
  V.store     = V[, 1];

  ptime = proc.time()[1]

  # Metropolis-within-Gibbs sampler

  for(j in 1:(N-1)){

    adj    = (j - 1) * thin;

    V.star = V.store; # proposal value

    aux    = sample(k1); # positions to be changed in the thining loop

    # Thining
    for (i in 1:thin) {

      iter = i + adj;

      if (iter %% printIter == 0) {
        cat(paste("Iteration", iter, ",", "Time elapsed",
                  round(as.numeric(proc.time()[1] - ptime) / 60, 2),
                  "minutes"), "\n")
      }

      f.store <- lpost(omega,
                       FZ,
                       k,
                       V.store,     # parameter
                       tau.store,   # parameter
                       tau.alpha,
                       tau.beta,
                       phi.store,   # parameter
                       phi.alpha,
                       phi.beta,
                       delta.store, # parameter
                       delta.alpha,
                       delta.beta,
                       P,
                       pdgrm,
                       degree,
                       db.list)

      ##############
      ### WEIGHT ###
      ##############

      #aux     = sample(k1);

      # tunning proposal distribution

      if(count < 0.30){ # increasing acceptance pbb

        sigma = sigma * 0.90; # decreasing proposal moves

      }else if(count > 0.50){ # decreasing acceptance pbb

        sigma = sigma * 1.1; # increasing proposal moves

      }

      count = 0; # ACCEPTANCE PROBABILITY

      for(g in 1:k1){

        pos         = aux[g];

        V.star[pos] = V.store[pos] + sigma * Zs[iter, g];

        f.V.star <- lpost(omega,
                          FZ,
                          k,
                          V.star, # proposal value
                          tau.store,
                          tau.alpha,
                          tau.beta,
                          phi.store,
                          phi.alpha,
                          phi.beta,
                          delta.store,
                          delta.alpha,
                          delta.beta,
                          P,
                          pdgrm,
                          degree,
                          db.list)

        # log posterior for previous iteration
        #f.V <- f.store;

        #Accept/reject

        alpha1 <- min(0, f.V.star - f.store); # log acceptance ratio

        if(Us[iter,g] < alpha1) {

          V.store[pos] <- V.star[pos];  # Accept W.star
          f.store      <- f.V.star;
          count        <- count + 1; # ACCEPTANCE PROBABILITY

        }else {

          V.star[pos] = V.store[pos]; # reseting proposal value

        }

      } # End updating weights

      count       = count / k1;
      Count[iter] = count; # Acceptance probability

      ###########
      ### phi ###
      ###########

      phi.store = stats::rgamma(1, shape = k/2 + phi.alpha,
                               rate = phi.beta * delta.store + t(V.store) %*% P %*% V.store / 2);

      #############
      ### delta ###
      #############

      delta.store = stats::rgamma(1, shape = phi.alpha + delta.alpha,
                                  rate = phi.beta * phi.store + delta.beta);

      ###########
      ### tau ###
      ###########

      # Directly sample tau from conjugate Inverse-Gamma density

      q.psd <- qpsd(omega, k, V.store, degree, db.list)$psd;
      m <- n - 2;
      q <- rep(NA, m);
      q[1] <- q.psd[1];
      q[m] <- q.psd[length(q.psd)];
      q[2 * 1:(m / 2 - 1)] <- q[2 * 1:(m / 2 - 1) + 1] <- q.psd[1:(m / 2 - 1) + 1];

      # Note the (n - 2) here - we remove the first and last terms
      tau.store <- 1 / stats::rgamma(1, tau.alpha + (n - 2) / 2,
                                      tau.beta + sum(pdgrm[2:(n - 1)] / q) / (2 * pi) / 2);

    }# End thining

  ######################
  ### Storing values ###
  ######################

  phi[j+1]   = phi.store;
  delta[j+1] = delta.store;
  tau[j+1]   = tau.store;
  V          = cbind(V, V.store);

  ### ###

}  # END: MCMC loop

  # Discarding burn-in period
  keep <- seq(round(burnin/thin) + 1, N);
  tau  <- tau[keep];
  phi  <- phi[keep];
  delta<- delta[keep];
  V    <- V[, keep];

  ##############################
  ### Compute log likelihood ###
  ##############################

  for(i in 1:length(keep)){

    ll.trace = c(ll.trace,
                 llike(omega, FZ, k, V[,i], tau[i], pdgrm, degree, db.list));
  }

  #fpsd.sample <- log.fpsd.sample <- matrix(NA, nrow = length(omega), ncol = length(keep));
  fpsd.sample <- log.fpsd.sample <- matrix(NA, nrow = length(omega) - 2, ncol = length(keep));

  # Store PSDs
  for (isample in 1:length(keep)) {
    q.psd <- qpsd(omega, k, V[, isample], degree, db.list); # db.list ADDED
    fpsd.sample[, isample] <- tau[isample] * q.psd$psd;
    #knots.trace[1:length(q.psd$knots), isample] <- q.psd$knots
    log.fpsd.sample[, isample] <- logfuller(fpsd.sample[, isample]); # Create transformed version
  }

  # Compute point estimates and 90% Pointwise CIs
  psd.median <- apply(fpsd.sample, 1, stats::median);
  psd.mean <- apply(fpsd.sample, 1, mean);
  psd.p05 <- apply(fpsd.sample, 1, stats::quantile, probs=0.05);
  psd.p95 <- apply(fpsd.sample, 1, stats::quantile, probs=0.95);

  # Transformed versions of these for uniform CI construction
  log.fpsd.s <- apply(log.fpsd.sample, 1, stats::median);
  log.fpsd.mad <- apply(log.fpsd.sample, 1, stats::mad);
  log.fpsd.help <- apply(log.fpsd.sample, 1, uniformmax);
  log.Cvalue <- stats::quantile(log.fpsd.help, 0.9);

  # Compute Uniform CIs
  psd.u95 <- exp(log.fpsd.s + log.Cvalue * log.fpsd.mad);
  psd.u05 <- exp(log.fpsd.s - log.Cvalue * log.fpsd.mad);

  ###########
  ### DIC ###
  ###########

  tau_mean = mean(tau);

  v_means = unname(apply(V, 1, mean));

  l = llike(omega, FZ, k, v = v_means,
            tau = tau_mean, pdgrm, degree, db.list);

  ls = apply(rbind(tau, V), 2, function(x){
             llike(omega, FZ, k, v = x[-1],
             tau = x[1], pdgrm, degree, db.list)});
  ls = unname(ls);

  # http://kylehardman.com/BlogPosts/View/6
  # DIC = -2 * (l - (2 * (l - mean(ls))));

  D_PostMean = -2 * l;
  D_bar      = -2 * mean(ls);
  pd         = D_bar - D_PostMean;

  DIC = list(DIC = 2 * D_bar - D_PostMean, pd = pd);

  # Compute periodogram
  N     = length(psd.median); # N = (n + 1) / 2 (ODD) or N = n / 2 + 1 (EVEN)
  pdgrm = (abs(stats::fft(data)) ^ 2 / (2 * pi * n))[1:N];

  # storing analysis specifications
  anSpecif = list(k = k, n = n, degree = degree, FZ = FZ, eqSpacedKnots = eqSpacedKnots);

  # List to output
  output = list(psd.median = psd.median * rescale ^ 2,
                psd.mean = psd.mean * rescale ^ 2,
                psd.p05 = psd.p05 * rescale ^ 2,
                psd.p95 = psd.p95 * rescale ^ 2,
                psd.u05 = psd.u05 * rescale ^ 2,
                psd.u95 = psd.u95 * rescale ^ 2,
                fpsd.sample = fpsd.sample * rescale ^ 2,
                anSpecif = anSpecif,
                n = n,
                tau = tau,
                phi = phi,
                delta = delta,
                V = V,
                ll.trace = ll.trace,
                pdgrm = pdgrm * rescale ^ 2,
                db.list = db.list,
                DIC = DIC,
                count = Count); # Acceptance probability

  class(output) = "psd"  # Assign S3 class to object

  return(output); # Return output

}  # Close function
