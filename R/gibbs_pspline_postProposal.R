#' This function takes into account a posterior pilot sample to calibrate
#'   the proposals for the weights
#'
#' @importFrom Rcpp evalCpp
#' @importFrom expm sqrtm
#' @importFrom stats median rnorm runif rgamma
#' @useDynLib psplinePsd, .registration = TRUE
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
                                       diffMatrixOrder = 3,
                                       psd,
                                       add = FALSE) {

  ### specific operations ###

  if(is.null(psd)){
    stop("Include an output from gibbs_psline function");
  }else{
    cat("Posterior samples used to calibrate the proposals for the weights", "\n")
  }

  k = psd$anSpecif$k;
  degree = psd$anSpecif$degree;

  cat(paste("Number of B-splines k=", k, sep=""), "\n");

  if( (Ntotal - burnin)/thin < k){
    stop("Change specifications: Either increase Ntotal or decrease thin
         parameters, otherwise the covariance matrix of the weights is singular")
  }

  ### ###

  if(!any(diffMatrixOrder == c(1,2,3))){
    stop("The order of the difference penalty matrix can only be 1, 2 or 3")
  }

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
  data    = data / rescale;  # Data now has standard deviation 1

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
  P       = diffMatrix(k-1, d = diffMatrixOrder); # Third order penalty
  P       = t(P) %*% P;
  epsilon = 1e-6; #1e-10;
  P       = P + epsilon * diag(dim(P)[2]); # P^(-1)=Sigma (Covariance matrix)

  ##################################
  ### Using gibbs_pspline output ###
  ##################################

  # Proposal calibration
  muV  = apply(psd$V, 1, mean); # mean vector - posterior samples
  covV = stats::cov(t(psd$V));  # covariance matrix - posterior samples
  sqrt.covV = expm::sqrtm(covV);

  # Starting values

  if( add ==TRUE){

    pilot_N = length(psd$tau);

    tau[1]   <- psd$tau[pilot_N];
    delta[1] <- psd$delta[pilot_N];
    phi[1]   <- psd$phi[pilot_N];
    v        <- psd$V[, pilot_N];

  }else{

    tau[1]   <- stats::median(psd$tau);  # stats::var(data) / (2 * pi);
    delta[1] <- stats::median(psd$delta);# delta.alpha / delta.beta;
    phi[1]   <- stats::median(psd$phi);  # phi.alpha/(phi.beta * delta[1]);
    v        <- apply(psd$V, 1, median);

  }

  v = solve(sqrt.covV) %*% (v - muV);

  V = matrix(v, ncol = 1);

  ###
  # Since the number of knots are fixed,
  #  the B-splines only need to be calculated once.
  newk     <- k - degree + 1;
  knots    <- seq(0,1, length = newk);
  db.list  <- dbspline(omega, knots, degree);
  ###

  # Store log likelihood
  ll.trace    <- rep(NA, Ntotal)
  ll.trace[1] <- llike(omega, FZ, k, sqrt.covV %*% V[, 1] + muV,
                       tau[1], pdgrm, degree, db.list)$llike

  Count    = NULL; # acceptance probability
  sigma    = 1;    # variance of proposal distb for weights
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

    v = sqrt.covV %*% V[, i] + muV;

    lp.list <- lpost(omega,
                     FZ,
                     k,
                     v,
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

    V.store = V[, i];

    V.star  = V.store;

    aux     = sample(k1);

    # tunning proposal distribution

    if(count < 0.30){ # increasing acceptance pbb

      sigma = sigma * 0.90; # decreasing proposal moves

    }else if(count > 0.50){ # decreasing acceptance pbb

      sigma = sigma * 1.1; # increasing proposal moves

    }
    #print(sigma)
    count   = 0; # acceptance probability

    for(j in 1:k1){
      #print(j)
      pos         = aux[j];

      #sigma       = Sigma[pos]; # Another approach, individual variance

      V.star[pos] = V.store[pos] + sigma * Zs[i,j];

      v = sqrt.covV %*% V.star + muV;

      lp.list <- lpost(omega,
                       FZ,
                       k,
                       v, # proposal value
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

      f.V.star <- lp.list$lp;

      # log posterior for previous iteration
      f.V <- f.store;

      #Accept/reject

      alpha1 <- min(0, f.V.star - f.V);  # log acceptance ratio

      if(Us[i,j] < alpha1) {

        V.store[pos] <- V.star[pos];  # Accept W.star
        f.store      <- f.V.star;
        count        <- count + 1; # ACCEPTANCE PROBABILITY

      }else {

        V.star[pos] = V.store[pos]; # reseting proposal value

      }

    } # End updating weights

    V        = cbind(V, V.store);
    count    = count / k1;
    Count[i] = count; # acceptance probability

    ###########
    ### phi ###
    ###########

    v        = sqrt.covV %*% V[, i+1] + muV;

    phi[i+1] = stats::rgamma(1, shape = k1/2 + phi.alpha,
                             rate = phi.beta * delta[i] + t(v) %*% P %*% v / 2);

    #############
    ### delta ###
    #############

    delta[i+1] = stats::rgamma(1, shape = phi.alpha + delta.alpha,
                               rate = phi.beta * phi[i+1] + delta.beta);

    ###########
    ### tau ###
    ###########

    # Directly sample tau from conjugate Inverse-Gamma density

    q.psd <- qpsd(omega, k, v, degree, db.list)$psd;
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

    ll.trace[i + 1] <- llike(omega, FZ, k, v, tau[i + 1], pdgrm,
                             degree, db.list)$llike;

  } # END: MCMC loop

  # Which iterations to keep
  keep  <- seq(burnin + 1, Ntotal, by = thin);
  tau   <- tau[keep];
  phi   <- phi[keep];
  delta <- delta[keep];
  V     <- V[, keep];
  ll.trace <- ll.trace[keep];

  V = apply(V,2, function(x) sqrt.covV %*% x + muV); # converting to actual V

  fpsd.sample <- log.fpsd.sample <- matrix(NA, nrow = length(omega) - 2, ncol = length(keep));
  # knots.trace <- matrix(NA, nrow = kmax, ncol = length(keep))

  # Store PSDs
  for (isample in 1:length(keep)) {
    #q.psd <- qpsd(omega, k, sqrt.covV %*% V[, isample] + muV, degree, db.list); # db.list ADDED
    q.psd <- qpsd(omega, k, V[, isample], degree, db.list); # db.list ADDED
    fpsd.sample[, isample] <- tau[isample] * q.psd$psd;
    #knots.trace[1:length(q.psd$knots), isample] <- q.psd$knots
    log.fpsd.sample[, isample] <- logfuller(fpsd.sample[, isample]); # Create transformed version
  }

  if(add == TRUE){

    cat("Pilot and current analyses have been merged", "\n");

    count = c(psd$count * k, Count); # k factor is cancelled return output

    # db.list is equal in both analyses

    delta = c(psd$delta, delta);

    aux   = psd$fpsd.sample / (rescale^2);

    fpsd.sample = cbind(aux, fpsd.sample);

    log.fpsd.sample = cbind(apply(aux, 2, logfuller), log.fpsd.sample);

    ll.trace = c(psd$ll.trace, ll.trace);

    # pdgrm is equal in both analyses

    phi = c(psd$phi, phi);

    # psd.mean, psd.median, psd.p05, psd.p95, psd.u95
    #   are calculated below from fpsd.sample

    tau = c(psd$tau, tau);

    V   = cbind(psd$V, V);

  }

  # Compute point estimates and 90% Pointwise CIs
  psd.median <- apply(fpsd.sample, 1, stats::median);
  psd.mean   <- apply(fpsd.sample, 1, mean);
  psd.p05    <- apply(fpsd.sample, 1, stats::quantile, probs=0.05);
  psd.p95    <- apply(fpsd.sample, 1, stats::quantile, probs=0.95);

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

  v_means = unname(apply(V, 1, mean)); # V was converted above

  l = llike(omega, FZ, k, v = v_means,
            tau = tau_mean, pdgrm, degree, db.list)$llike;

  ls = apply(rbind(tau, V), 2, function(x){
             llike(omega, FZ, k, v = x[-1],
             tau = x[1], pdgrm, degree, db.list)$llike});

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
  anSpecif = list(k = k, n = n, degree = degree, FZ = FZ);

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
                count = Count/k); # Acceptance probability

  class(output) = "psd"  # Assign S3 class to object

  return(output); # Return output

  }  # Close function
