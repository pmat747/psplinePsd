#' This function takes into account a posterior pilot sample to calibrate
#'   the proposals for the weights
#'
#' @importFrom Rcpp evalCpp
#' @keywords internal
gibbs_pspline_postProposal <- function(data,
                                       Ntotal,
                                       burnin,
                                       thin = 1,
                                       tau.alpha = 0.001,
                                       tau.beta = 0.001,
                                       phi.alpha = 1,
                                       phi.beta = 1e-5,
                                       delta.alpha,
                                       delta.beta,
                                       k = NULL,
                                       degree = 3,
                                       psd) {

  ### specific operations ###

  if(is.null(psd)){
    stop("include an output from gibbs_psline function");
  }

  if( (Ntotal - burnin)/thin < k){
    stop("Change specifications: Either increase Ntotal or decrease thin
         parameters, otherwise the covariance matrix of the weights is singular")
  }

  k = psd$k;

  ### ###

  n <- length(data);

  if (n %% 2 != 0) stop("this version of bsplinePsd must have n even")

  # Tolerance for mean centering
  tol <- 1e-4;

  # Mean center
  if (abs(mean(data)) > tol) {
    data <- data - mean(data)
    warning("Data has been mean-centred")
  }

  # Optimal rescaling to prevent numerical issues
  rescale = stats::sd(data);
  data = data / rescale;  # Data now has standard deviation 1

  if (burnin >= Ntotal) stop("burnin must be less than Ntotal")
  if (any(c(Ntotal, burnin, thin) %% 1 != 0 || any(c(Ntotal, burnin, thin) < 0)))
  if (any(c(tau.alpha, tau.beta) <= 0)) stop(" tau.alpha and tau.beta must be strictly positive")
  if (any(c(Ntotal, thin) %% 1 != 0) || any(c(Ntotal, thin) <= 0)) stop("Ntotal must be strictly positive integers")
  if ((burnin %% 1 != 0) || (burnin < 0)) stop("burnin must be a non-negative integer")

  FZ <- fast_ft(data); # FFT data to frequency domain.  NOTE: Must be mean-centred.

  pdgrm <- abs(FZ) ^ 2; # Periodogram: NOTE: the length is n here.

  omega <- 2 * (1:(n / 2 + 1) - 1) / n; # Frequencies on unit interval
  lambda <- pi * omega; # Angular frequencies on [0, pi]

  # Open objects for storage
  tau   <- rep(NA, Ntotal);
  phi   <- rep(NA, Ntotal); # scalar in Multivariate Normal
  delta <- rep(NA, Ntotal);

  # Difference Matrix
  P       = diffMatrix(k-1, d = 3); # Third order penalty
  P       = t(P) %*% P;
  epsilon = 1e-6; #1e-10;
  P       = P + epsilon * diag(dim(P)[2]); # P^(-1)=Sigma (Covariance matrix)

  ##################################
  ### Using gibbs_pspline output ###
  ##################################

  # Starting values
  tau[1]   <- stats::median(psd$tau);  # stats::var(data) / (2 * pi);
  delta[1] <- stats::median(psd$delta);# delta.alpha / delta.beta;
  phi[1]   <- stats::median(psd$phi);  # phi.alpha/(phi.beta * delta[1]);

  # starting value for the weigths
  muW  = apply(psd$W, 1, mean); # mean vector - posterior samples
  covW = stats::cov(t(psd$W));  # covariance matrix - posterior samples
  sqrt.covW = expm::sqrtm(covW);
  #Sigma = sqrt(diag(covW));

  w = pdgrm / sum(pdgrm);
  w = w[round(seq(1, length(w), length = k))];
  w[which(w==0)] = 1e-50; # prevents errors when there are zeros
  w = w/sum(w);
  w = w[-k];
  w = log(w / (1 - sum(w)));

  w = (w - muW) %*% solve(sqrt.covW);

  W = matrix(w, ncol = 1);

  ###
  # Since the number of knots are fixed,
  #  the B-splines only need to be calculated once.
  newk     <- k - degree + 1;
  knots    <- seq(0,1, length = newk);
  db.list  <- dbspline(omega, knots, degree);
  ###

  # Store log likelihood
  ll.trace    <- rep(NA, Ntotal)
  ll.trace[1] <- llike(omega, FZ, k, sqrt.covW %*% W[, 1] + muW,
                       tau[1], pdgrm, degree, db.list)$llike

  Count    = NULL; # acceptance probability
  sigma    = 1;    # variance of proposal distb for weigths
  count    = 0.4;  # starting value for acc pbb - optimal value

  k1 = k - 1;

  # Random values
  Zs = stats::rnorm((Ntotal-1)*k1);
  Zs = matrix(Zs, nrow = Ntotal-1, ncol = k1);
  Us = log(stats::runif((Ntotal-1)*k1, 0, 1));
  Us = matrix(Us, nrow = Ntotal-1, ncol = k1);

  ptime = proc.time()[1]

  # Metropolis-within-Gibbs sampler
  for (i in 1:(Ntotal-1)) {

    if (i %% 100 == 0) {
      cat(paste("Iteration", i, ",", "Time elapsed",
                round(as.numeric(proc.time()[1] - ptime) / 60, 2),
                "minutes"), "\n")
    }

    w = sqrt.covW %*% W[, i] + muW;

    lp.list <- lpost(omega,
                     FZ,
                     k,
                     w,
                     tau[i],
                     tau.alpha,
                     tau.beta,
                     phi[i],
                     phi.alpha,
                     phi.beta,
                     delta[i],
                     delta.alpha,
                     delta.beta,
                     P,
                     pdgrm,
                     degree,
                     db.list)

    f.store <- lp.list$lp;

    ##############
    ### WEIGHT ###
    ##############

    W.store = W[, i];

    W.star  = W.store;

    aux     = sample(k1);

    # tunning proposal distribution

    if(count < 0.30){ # increasing acceptance pbb

      sigma = sigma * 0.90; # decreasing proposal moves

    }else if(count > 0.50){ # decreasing acceptance pbb

      sigma = sigma * 1.1; # increasing proposal moves

    }
    #print(sigma)
    count   = 0;# ACCEPTANCE PROBABILITY

    for(j in 1:k1){
      #print(j)
      pos         = aux[j];

      #sigma       = Sigma[pos]; # Another approach, individual variance

      W.star[pos] = W.store[pos] + sigma * Zs[i,j];

      w = sqrt.covW %*% W.star + muW;

      lp.list <- lpost(omega,
                       FZ,
                       k,
                       w, # proposal value
                       tau[i],
                       tau.alpha,
                       tau.beta,
                       phi[i],
                       phi.alpha,
                       phi.beta,
                       delta[i],
                       delta.alpha,
                       delta.beta,
                       P,
                       pdgrm,
                       degree,
                       db.list)

      f.W.star <- lp.list$lp;

      # log posterior for previous iteration
      f.W <- f.store;

      #Accept/reject

      alpha1 <- min(0, f.W.star - f.W);  # log acceptance ratio

      if(Us[i,j] < alpha1) {

        W.store[pos] <- W.star[pos];  # Accept W.star
        f.store      <- f.W.star;
        count        <- count + 1; # ACCEPTANCE PROBABILITY

      }else {

        W.star[pos] = W.store[pos]; # reseting proposal value

      }

    } # End updating weights

    W        = cbind(W, W.store);
    count    = count / k1;
    Count[i] = count; # acceptance probability

    ###########
    ### phi ###
    ###########

    w        = sqrt.covW %*% W[, i+1] + muW;

    phi[i+1] = stats::rgamma(1, shape = k1/2 + phi.alpha,
                             rate = phi.beta * delta[i] + t(w) %*% P %*% w / 2);

    #############
    ### delta ###
    #############

    delta[i+1] = stats::rgamma(1, shape = phi.alpha + delta.alpha,
                               rate = phi.beta * phi[i+1] + delta.beta);

    ###########
    ### tau ###
    ###########

    # Directly sample tau from conjugate Inverse-Gamma density

    q.psd <- qpsd(omega, k, w, degree, db.list)$psd;
    m     <- n - 2;
    q     <- rep(NA, m);
    q[1]  <- q.psd[1];
    q[m]  <- q.psd[length(q.psd)];
    q[2 * 1:(m / 2 - 1)] <- q[2 * 1:(m / 2 - 1) + 1] <- q.psd[1:(m / 2 - 1) + 1];

    # Note the (n - 2) here - we remove the first and last terms
    tau[i+1] <- 1 / stats::rgamma(1, tau.alpha + (n - 2) / 2,
                                  tau.beta + sum(pdgrm[2:(n - 1)] / q) / (2 * pi) / 2);

    ##############################
    ### Compute log likelihood ###
    ##############################

    ll.trace[i + 1] <- llike(omega, FZ, k, w, tau[i + 1], pdgrm,
                             degree, db.list)$llike;

  } # END: MCMC loop

  # Which iterations to keep
  keep  <- seq(burnin + 1, Ntotal, by = thin);
  tau   <- tau[keep];
  phi   <- phi[keep];
  delta <- delta[keep];
  W     <- W[, keep];
  ll.trace <- ll.trace[keep];

  fpsd.sample <- log.fpsd.sample <- matrix(NA, nrow = length(omega) - 2, ncol = length(keep));
  # knots.trace <- matrix(NA, nrow = kmax, ncol = length(keep))

  # Store PSDs
  for (isample in 1:length(keep)) {
    q.psd <- qpsd(omega, k, sqrt.covW %*% W[, isample] + muW, degree, db.list); # db.list ADDED
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

  # Compute periodogram
  N     = length(psd.median); # N = (n + 1) / 2 (ODD) or N = n / 2 + 1 (EVEN)
  pdgrm = (abs(stats::fft(data)) ^ 2 / (2 * pi * n))[1:N];

  ###########
  ### DIC ###
  ###########

  tau_mean = mean(tau);

  w_means = unname(apply(W, 1, mean));

  l = llike(omega, FZ, k, w = sqrt.covW %*% w_means + muW,
            tau = tau_mean, pdgrm, degree, db.list)$llike;

  ls = apply(rbind(tau, W), 2, function(x){
    llike(omega, FZ, k, w = sqrt.covW %*% x[-1] + muW,
          tau = x[1], pdgrm, degree, db.list)$llike});
  ls = unname(ls);

  # http://kylehardman.com/BlogPosts/View/6
  # DIC = -2 * (l - (2 * (l - mean(ls))));

  D_PostMean = -2 * l;
  D_bar      = -2 * mean(ls);
  pd         = D_bar - D_PostMean;

  DIC = list(DIC = 2 * D_bar - D_PostMean, pd = pd);

  # List to output
  output = list(psd.median = psd.median * rescale ^ 2,
                psd.mean = psd.mean * rescale ^ 2,
                psd.p05 = psd.p05 * rescale ^ 2,
                psd.p95 = psd.p95 * rescale ^ 2,
                psd.u05 = psd.u05 * rescale ^ 2,
                psd.u95 = psd.u95 * rescale ^ 2,
                fpsd.sample = fpsd.sample,
                k = k,
                tau = tau,
                phi = phi,
                delta = delta,
                W = W,
                ll.trace = ll.trace,
                pdgrm = pdgrm * rescale ^ 2,
                n = n,
                db.list = db.list,
                DIC = DIC,
                count = Count/k); # Acceptance probability

  class(output) = "psd"  # Assign S3 class to object

  return(output); # Return output

  }  # Close function
