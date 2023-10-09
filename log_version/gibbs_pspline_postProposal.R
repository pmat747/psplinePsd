#' This function takes into account a posterior pilot sample to calibrate
#'   the proposals for the weights
#'
#' @importFrom Rcpp evalCpp
#' @importFrom expm sqrtm
#' @importFrom stats median rnorm runif rgamma
#' @importFrom fda create.bspline.basis bsplinepen
#' @useDynLib psplinePsd, .registration = TRUE
#' @keywords internal


# Adaptive Metropolis Hasting

gibbs_pspline_postProposal <- function(data,
                                       Ntotal,
                                       burnin,
                                       thin,
                                       tau.alpha,
                                       tau.beta,
                                       phi.alpha,
                                       phi.beta,
                                       delta.alpha,
                                       delta.beta,
                                       k,
                                       eqSpacedKnots,
                                       nfreqbin,
                                       wfreqbin,
                                       degree,
                                       diffMatrixOrder,
                                       printIter,
                                       amh) {
  
  if(burnin >= Ntotal) stop("burnin must be less than Ntotal");
  if(any(c(Ntotal, thin) %% 1 != 0) || any(c(Ntotal, thin) <= 0)) stop("Ntotal and thin must be strictly positive integers");
  if((burnin %% 1 != 0) || (burnin < 0)) stop("burnin must be a non-negative integer");
  if(any(c(tau.alpha, tau.beta) <= 0)) stop("tau.alpha and tau.beta must be strictly positive");
  if(any(c(phi.alpha, phi.beta) <= 0)) stop("phi.alpha and phi.beta must be strictly positive");
  if(any(c(delta.alpha, delta.beta) <= 0)) stop("delta.alpha and delta.beta must be strictly positive");
  
  if(!is.logical(eqSpacedKnots))stop("eqSpacedKnots must be a logical value");
  
  n <- length(data);
  
  # Which boundary frequencies to remove from likelihood computation and tau sample
  if (n %% 2) { # Odd length time series
    bFreq <- 1; # Remove first
  }
  else {  # Even length time series
    bFreq <- c(1, n); # Remove first and last
  }

  if(is.null(k)){
    
    k = min(round(n/4), 40);
    cat(paste("The number of B-splines has not been specified. Therefore, k=min(round(n/4), 40)=", k, sep=""), "\n");
    
  }else{
    
    if((k %% 1 !=0) || (k <= 0))stop("k must be strictly positive integer");
    if((k < degree + 2 ))stop("k must be at least degree + 2");
    
  }

  if((Ntotal - burnin)/thin < k)stop("Change corresponding specifications in order to have (Ntotal-burnin)/thin > k");
  
  if((printIter<=0) || (printIter %% 1 != 0) )stop("printIter must be a positive integer value");
  
  if(!(k-2 >= diffMatrixOrder))stop("diffMatrixOrder must be lower than or equal to k-2");
  
  # Tolerance for mean centering
  tol <- 1e-4;
  
  # Mean center
  if (abs(mean(data)) > tol) {
    data <- data - mean(data);
    warning("Data has been mean-centred");
  }
  
  # Optimal rescaling to prevent numerical issues
  rescale = stats::sd(data);
  data    = data / rescale;  # Data now has standard deviation 1
  
  FZ     <- fast_ft(data); # FFT data to frequency domain.  NOTE: Must be mean-centred.
  
  pdgrm  <- abs(FZ) ^ 2; # Periodogram: NOTE: the length is n here.
  
  omega  <- 2 * (1:(n / 2 + 1) - 1) / n;  # Frequencies on unit interval
  lambda <- pi * omega; # Angular frequencies on [0, pi]
  
  # log-psd AR(1)
  orderAR = 1;
  spec_ar = stats::spectrum(data, n.freq = length(FZ), method = "ar", 
                            plot = FALSE, order = orderAR)$spec;
  
  spec_ar = log(spec_ar);
  
  N = round(Ntotal/thin); # number of samples to be stored
  
  # Open objects for storage
  phi   <- rep(NA, N); # scalar in Multivariate Normal
  delta <- rep(NA, N); # Hyperparameter
  tau   <- rep(NA, N); # normalising constant
  
  # Starting values
  delta[1] <- delta.alpha/delta.beta;
  phi[1]   <- phi.alpha/(phi.beta * delta[1]);
  tau[1]   <- stats::var(data) / (2 * pi);
  
  # Fixed knot number & location => a single calculation of B-spline densities
  knots = knotLoc(data = data, k = k, degree = degree, eqSpaced = eqSpacedKnots,
                  nfreqbin = nfreqbin, wfreqbin = wfreqbin); # length(knots) = k - degree + 1
 # k     = length(knots) + degree - 1; # Correcting number of B-splines

  db.list <- dbspline(omega, knots, degree);
  db.list <- apply(db.list, 2, function(x) (x - mean(x))/sd(x)); # Standarization
  
  # starting value for the weights
  w = pdgrm / sum(pdgrm); # poor alternative: w=rep(0,k-1);
  w = w[round(seq(1, length(w), length = k))];
  w[which(w==0)] = 1e-50; # prevents errors when there are zeros
  v = w/sum(w); # log(v)?
  V = matrix(v, ncol = 1);
  
  # Difference Matrix
  if(eqSpacedKnots == TRUE){
    
    P = diffMatrix(k, d = diffMatrixOrder);
    P = t(P) %*% P;
    
  }else{
    
    if(degree <= diffMatrixOrder){
      stop("bsplinepen function: penalty order must be lower than the bspline density degree");
    }
    
    basisobj = fda::create.bspline.basis(c(0, 1), norder = degree + 1,
                                        nbasis = k , breaks = knots);
    
    P = fda::bsplinepen(basisobj, Lfdobj = diffMatrixOrder);
    
    P = P/norm(P);
  }
  
  epsilon = 1e-6; #1e-6;
  P       = P + epsilon * diag(dim(P)[2]); # P^(-1)=Sigma (Covariance matrix)

  ll.trace   = NULL;  # Store log likelihood
  Count      = NULL;  # acceptance probability
  sigta      = 1;     # variance of proposal distb for tau
  sigtaStore = NULL;
  count      = NULL;  # 
  count_tau  = 0.4;   # counting for tau
  count_tau_alpha = 0;# counting for tau
  Count_tau = NULL;
  Count_tau_alpha = NULL;
  
  # Random values
  
  # V (weights)
  Uv_am = stats::runif((Ntotal-1)*thin, 0, 1); # AMH step
  Uv_am = matrix(Uv_am, nrow = Ntotal-1, ncol = thin);
  Uv = log(stats::runif((Ntotal-1)*thin, 0, 1));  # Acceptance step
  Uv = matrix(Uv, nrow = Ntotal-1, ncol = thin);
  
  # tau
  Zt = stats::rnorm((Ntotal-1)*thin);            # Proposals
  Zt = matrix(Zt, nrow = Ntotal-1, ncol = thin);
  Ut = log(stats::runif((Ntotal-1)*thin, 0, 1)); # Acceptance step
  Ut = matrix(Ut, nrow = Ntotal-1, ncol = thin);
  
  # Initial values for the proposals
  phi.store   = phi[1];
  tau.store   = tau[1];
  delta.store = delta[1];
  V.store     = V[, 1];
  
  covObj = updateCov(X = V.store, covObj = NA);
  Ik     = 0.1^2 * diag(rep(1, k)/k); # matrix used in AMH step
  c_amh  = 2.38^2 / k; # constant used in AMH step
  
  # Time
  ptime = proc.time()[1];
  
  # Calculating likelihood for starting values
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
                   db.list,
                   spec_ar);
  
  ll.trace = c(ll.trace, f.store$ll); # Storing likelihood
  
  # Metropolis-within-Gibbs sampler
  
  for(j in 1:(N-1)){
    
    adj    = (j - 1) * thin;
    
    # Thining
    for(i in 1:thin) {
      
      iter = i + adj;
      
      if(iter %% printIter == 0) {
        cat(paste("Iteration", iter, ",", "Time elapsed",
                  round(as.numeric(proc.time()[1] - ptime) / 60, 2),
                  "minutes"), "\n");
      }
      
      ##############
      ### WEIGHT ###
      ##############
      
      if( j <= 2*k ){
        
        V.star = mvtnorm::rmvnorm(1, mean = V.store, sigma = Ik);
        
      }else{
        
        if(Uv_am[iter, i] < 0.05){
          
          V.star = mvtnorm::rmvnorm(1, mean = V.store, sigma = Ik);
          
        }else{
          
          V.star = mvtnorm::rmvnorm(1, mean = V.store, sigma = c_amh * S);
          
        }
      }
      
      V.star = t(V.star);
      
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
                        db.list,
                        spec_ar);
      
      #Accept/reject
      alpha1 <- min(0, f.V.star$lp - f.store$lp); # log acceptance ratio
      
      if(Uv[iter, i] < alpha1){
        
        V.store = V.star; # Accept W.star
        f.store = f.V.star;
        count   = c(count, 1); # acceptance probability
        
      }else{
        
        count  = c(count, 0);
        
      }
      
      ###########
      ### tau ###
      ###########
      
      tau.star = tau.store + sigta * Zt[iter, i];
      
      f.tau.star <- lpost(omega,
                          FZ,
                          k,
                          V.store, 
                          tau.star, # proposal value
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
                          db.list,
                          spec_ar);
      
      alpha_tau <- min(0, f.tau.star$lp - f.store$lp); # log acceptance ratio
      
      if(Ut[iter, i] < alpha_tau){
        
        tau.store <- tau.star;      # Accept tau.star
        f.store   <- f.tau.star;
        count_tau <- count_tau + 1; # acceptance probability
        
      }
      
    } # Thinning # End updating weights
    
    Count[iter] = mean(count); # Acceptance probability
    count = NULL; # reseting count
    
    count_tau = count_tau / k;
    Count_tau = c(Count_tau, count_tau);
    
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
    
    ######################
    ### Storing values ###
    ######################
    
    phi[j+1]   = phi.store;
    delta[j+1] = delta.store;
    tau[j+1]   = tau.store;
    V          = cbind(V, V.store);
    
    # Calculating likelihoods of the new values
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
                     db.list,
                     spec_ar);
    
    ll.trace   = c(ll.trace, f.store$ll); 
    
    # Updating covariance matrix
    covObj = updateCov(X = V.store, covObj = covObj);
    S      = covObj$cov;
    
  }  # END: MCMC loop
  
  print("Processing posterior samples")
  
  # Discarding burn-in period
  keep     <- seq(round(burnin/thin) + 1, N);
  tau      <- tau[keep];
  phi      <- phi[keep];
  delta    <- delta[keep];
  V        <- V[, keep];
  ll.trace <- ll.trace[keep]; 
  
  ####################
  ### Post process ###
  ####################
  
  log.fpsd.sample <- matrix(NA, nrow = length(omega), ncol = length(keep));
  
  ptime = proc.time()[1];
  
  spec_ar = stats::spectrum(data, n.freq = length(omega), method = "ar", 
                            plot = FALSE, order = orderAR)$spec;
  
  spec_ar = log(spec_ar);
  
  # Store PSDs
  for (isample in 1:length(keep)){
    
    q.psd <- qpsd(omega, k, V[, isample], degree, db.list);
    log.fpsd.sample[, isample] <- tau[isample] + q.psd + spec_ar; 
   
  }
  
  cat(paste("qpsd ", 
            round(as.numeric(proc.time()[1] - ptime) / 60, 2),
            "minutes"), "\n");
  
  # Compute point estimates and 90% Pointwise CIs
  #  psd.median <- apply(fpsd.sample, 1, stats::median);
  #  psd.mean <- apply(fpsd.sample, 1, mean);
  #  psd.p05 <- apply(fpsd.sample, 1, stats::quantile, probs=0.05);
  #  psd.p95 <- apply(fpsd.sample, 1, stats::quantile, probs=0.95);
  
  # Transformed versions of these for uniform CI construction
  #  log.fpsd.sample = fpsd.sample;
  ptime = proc.time()[1]
  log.psd.median <- apply(log.fpsd.sample, 1, stats::median);
  cat(paste("Median ", 
            round(as.numeric(proc.time()[1] - ptime) / 60, 4),
            "minutes"), "\n");
  ptime = proc.time()[1]
  log.psd.mad    <- apply(log.fpsd.sample, 1, stats::mad);
  cat(paste("Mad ", 
            round(as.numeric(proc.time()[1] - ptime) / 60, 4),
            "minutes"), "\n");
  ptime = proc.time()[1]
  log.psd.help   <- apply(log.fpsd.sample, 1, uniformmax);
  cat(paste("uniformmax ", 
            round(as.numeric(proc.time()[1] - ptime) / 60, 4),
            "minutes"), "\n");
  ptime = proc.time()[1]
  log.Cvalue      <- stats::quantile(log.psd.help, 0.9);
  cat(paste("uniformax ", 
            round(as.numeric(proc.time()[1] - ptime) / 60, 4),
            "minutes"), "\n");
  
  # Compute Uniform CIs
  ptime = proc.time()[1]
  psd.u95 <- log.psd.median + log.Cvalue * log.psd.mad;
  cat(paste("u95 ", 
            round(as.numeric(proc.time()[1] - ptime) / 60, 4),
            "minutes"), "\n");
  ptime = proc.time()[1]
  psd.u05 <- log.psd.median - log.Cvalue * log.psd.mad;
  cat(paste("u05 ", 
            round(as.numeric(proc.time()[1] - ptime) / 60, 4),
            "minutes"), "\n");
  
  ###########
  ### DIC ###
  ###########
  #  ptime = proc.time()[1]
  #  tau_mean = mean(tau);
  
  #  v_means = unname(apply(V, 1, mean));
  
  #  l = llike(omega, FZ, k, v = v_means,
  #            tau = tau_mean,tau.alpha, pdgrm, degree, db.list);
  
  #  ls = apply(rbind(tau, V), 2, function(x){
  #    llike(omega, FZ, k, v = x[-1],
  #          tau = x[1], tau.alpha, pdgrm, degree, db.list)});
  #  ls = unname(ls);
  
  # http://kylehardman.com/BlogPosts/View/6
  # DIC = -2 * (l - (2 * (l - mean(ls))));
  
  #  D_PostMean = -2 * l;
  #  D_bar      = -2 * mean(ls);
  #  pd         = D_bar - D_PostMean;
  
  #  DIC = list(DIC = 2 * D_bar - D_PostMean, pd = pd);
  #  cat(paste("DIC ", 
  #            round(as.numeric(proc.time()[1] - ptime) / 60, 4),
  #            "minutes"), "\n");
  
  # Compute periodogram
  N     = length(log.psd.median); # N = (n + 1) / 2 (ODD) or N = n / 2 + 1 (EVEN)
  pdgrm = (abs(stats::fft(data)) ^ 2 / (2 * pi * n))[1:N];
  
  # storing analysis specifications
  anSpecif = list(k = k, n = n, degree = degree, FZ = FZ,
                  eqSpacedKnots = eqSpacedKnots, diffMatrixOrder = diffMatrixOrder,
                  tau.alpha = tau.alpha, tau.beta = tau.beta,
                  phi.alpha = phi.alpha, phi.beta = phi.beta,
                  delta.alpha = delta.alpha, delta.beta = delta.beta,
                  S = S);
  
  # List to output
  output = list(#psd.median = psd.median * rescale ^ 2,
    #psd.mean = psd.mean * rescale ^ 2,
    #psd.p05 = psd.p05 * rescale ^ 2,
    #psd.p95 = psd.p95 * rescale ^ 2,
    log.psd.median = log.psd.median + 2*log(rescale),#fpsd.sample * rescale ^ 2,
    psd.u05 = psd.u05 + 2*log(rescale), # psd.u05 * rescale ^ 2,
    psd.u95 = psd.u95 + 2*log(rescale), # psd.u95 * rescale ^ 2
    fpsd.sample = log.fpsd.sample + 2*log(rescale),#fpsd.sample * rescale ^ 2,
    anSpecif = anSpecif,
    n = n,
    tau = tau,
    phi = phi,
    delta = delta,
    V = V,
    ll.trace = ll.trace,
    pdgrm = pdgrm * rescale ^ 2,
    db.list = db.list,
    # DIC = DIC,
    count = Count,
    Count_tau = Count_tau, # Acceptance probability
    sigtaStore = sigtaStore);
  
  class(output) = "psd"  # Assign S3 class to object
  print("sdjsajdkjsafdkjsa")
  return(output); # Return output
  
}  # Close function

#https://github.com/santafemha/samMCMC/blob/master/R/updateCov.
updateCov <- function(X,covObj=NA) {
  
  if(all(is.na(covObj))) {
    covObj <- list(mean=X,cov=NA,n=1)
    return(covObj)
  }
  
  covObj$n <- covObj$n + 1 # Update number of observations
  
  if(covObj$n==2) {
    X1 <- covObj$mean
    covObj$mean <- X1/2 + X/2
    dX1 <- X1 - covObj$mean
    dX2 <- X - covObj$mean
    covObj$cov <- tcrossprod(dX1,dX1) + tcrossprod(dX2,dX2)
    return(covObj)
  }
  
  dx <- covObj$mean - X # previous mean minus new X
  covObj$cov <- covObj$cov * (covObj$n-2)/(covObj$n-1) + tcrossprod(dx,dx)/covObj$n
  covObj$mean <- covObj$mean*(covObj$n-1)/covObj$n + X/covObj$n
  return(covObj)
}