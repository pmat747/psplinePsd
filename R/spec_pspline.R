#' @title Spectrum estimate via p-splines
#' @description This function allows to estimate the spectrum
#' @param data numeric vector
#' @param l length of the data subset to be analysed
#' @param p overlapping percentage of the data subsets
#' @param eq if is \code{TRUE}, last data subset has length \code{l}, even though the percentage is not \code{p}
#' @param Ntotal1 total number of iterations to run the pilot Markov chain
#' @param burnin1 number of initial iterations to be discarded in the pilot mcmc run
#' @param thin1 thinning number in the pilot mcmc run (post-processing)
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (post-processing)
#' @param tau.alpha,tau.beta prior parameters for tau (Inverse-Gamma)
#' @param phi.alpha,phi.beta prior parameters for phi (Gamma)
#' @param delta.alpha,delta.beta prior parameters for delta (Gamma)
#' @param k number of B-spline densities in the mixture
#' @param eqSpacedKnots logical value indicating whether the knots are equally spaced or defined according to the periodogram
#' @param degree positive integer specifying the degree of the B-spline densities (default is 3)
#' @param diffMatrixOrder positive integer specifying the order of the difference penalty matrix in the P-splines (default is 2)
#' @param printIter positive integer specifying the periodicity of the iteration number to be printed on screen (default 100)
#' @param recycl if \code{TRUE} the mcmc analysis for a data subset is used to calibrate the proposal for the succesor data subset.  It works for large \code{p} values.
#' @param likePlot if \code{TRUE}, a likelihood traceplot is generated for each data subset
#' @return A list with S3 class `psds' containing the the power spectral estimates for the data subsets and and a list with relevant information about the analysis
#' @seealso \link{gibbs_pspline}
#' @examples
#' \dontrun{
#'
#' set.seed(1)
#'
#' # Generate AR(1) data with rho = 0.9
#' n = 128
#' data = arima.sim(n, model = list(ar = 0.9));
#' data = data - mean(data);
#'
#' # Spectrum estimate via p-splines (may take some time)
#' spec = spec_pspline(data, l=50, p=90, Ntotal1=2000, burnin1=500, thin1=5,
#'                     Ntotal=1000, burnin=500, thin=5, k=30, recycl = TRUE);
#'
#' image(spec) # Plot log PSD (see documentation of image.plot.psd)
#'
#' }
#' @importFrom Rcpp evalCpp
#' @useDynLib psplinePsd, .registration = TRUE
#' @export
spec_pspline = function(data,
                        l,
                        p,
                        eq = TRUE,
                        Ntotal1,
                        burnin1,
                        thin1,
                        Ntotal,
                        burnin,
                        thin,
                        tau.alpha = 0.001,
                        tau.beta = 0.001,
                        phi.alpha = 1,
                        phi.beta = 1,
                        delta.alpha = 1e-04,
                        delta.beta = 1e-04,
                        k = NULL,
                        eqSpacedKnots = FALSE,
                        degree = 3,
                        diffMatrixOrder = 3,
                        printIter = 1000,
                        recycl = FALSE,
                        likePlot = FALSE){

  if (l %% 2 != 0) stop("this version of bsplinePsd must have l even")

  out   = list();

  index = ints(n = length(data), l = l, p = p, eq = eq);

  N     = dim(index)[1];

  if(is.null(k) && (eq == FALSE)){

    # if k has not been specified & the last data subset has a
    #  different length (eq=FALSE), so the las analysis might be
    #  performed with a different number of p-splines
    #  which might cause problem if recycl is TRUE

    sn = index[2,1] - index[N,1] + 1;

    if( min(c(sn/4, 40)) != min(c(l/4,40)) ){
      k = round(min(c(l/4,40)));
      warning(paste("k has been set to ", k, ", thus the last subset analysis uses the same number of B-splines as the rest", sep = ""),"\n");
    }
  }

  ptime = proc.time()[1];

  if(recycl == FALSE){

    for(i in 1:N){

      cat("Processing Interval number ", i, "\n");

      auxData = data[index[i,1]:index[i,2]];
      auxData = auxData - mean(auxData); # Mean center

      pilotmcmc = gibbs_pspline(data = auxData, Ntotal = Ntotal1, burnin = burnin1,
                                thin = thin1,
                                tau.alpha = tau.alpha, tau.beta = tau.beta,
                                phi.alpha = phi.alpha, phi.beta = phi.beta,
                                delta.alpha = delta.alpha, delta.beta = delta.beta,
                                k = k, eqSpacedKnots = eqSpacedKnots, degree = degree, diffMatrixOrder = diffMatrixOrder,
                                printIter = printIter, psd = NULL, add = FALSE);

      mcmc = gibbs_pspline(data = auxData, Ntotal = Ntotal, burnin = burnin,
                           thin = thin,
                           tau.alpha = tau.alpha, tau.beta = tau.beta,
                           phi.alpha = phi.alpha, phi.beta = phi.beta,
                           delta.alpha = delta.alpha, delta.beta = delta.beta,
                           k = k, eqSpacedKnots = eqSpacedKnots,
                           degree = degree, diffMatrixOrder = diffMatrixOrder,
                           printIter = printIter,
                           psd = pilotmcmc, add = FALSE) ;

      out[[i]] = mcmc$fpsd.sample;

      if(likePlot == TRUE){
        graphics::par();
        stats::ts.plot(mcmc$ll.trace, xlab = "Iteration", ylab = "log-likelihood",
                       main = paste("subset ", i, sep=""));
      }
    }

  }else{

    cat("Each mcmc analysis is used for the next analysis", "\n")

    auxData = data[index[1,1]:index[1,2]];
    auxData = auxData - mean(auxData); # Mean center

    cat("Initial pilot posterior sampling", "\n")
    mcmc = gibbs_pspline(data = auxData, Ntotal = Ntotal1, burnin = burnin1,
                         thin = thin1,
                         tau.alpha = tau.alpha, tau.beta = tau.beta,
                         phi.alpha = phi.alpha, phi.beta = phi.beta,
                         delta.alpha = delta.alpha, delta.beta = delta.beta,
                         k = k, eqSpacedKnots = eqSpacedKnots,
                         degree = degree, diffMatrixOrder = diffMatrixOrder,
                         printIter = printIter, psd = NULL, add = FALSE);

    for(i in 1:N){

      cat("Processing Interval number ", i, "\n");

      auxData = data[index[i,1]:index[i,2]];
      auxData = auxData - mean(auxData); # Mean center

      mcmc = gibbs_pspline(data = auxData, Ntotal = Ntotal, burnin = burnin,
                           thin = thin,
                           tau.alpha = tau.alpha, tau.beta = tau.beta,
                           phi.alpha = phi.alpha, phi.beta = phi.beta,
                           delta.alpha = delta.alpha, delta.beta = delta.beta,
                           k = k, eqSpacedKnots = eqSpacedKnots,
                           degree = degree, diffMatrixOrder = diffMatrixOrder,
                           printIter = printIter, psd = mcmc, add = FALSE);

      out[[i]] = mcmc$fpsd.sample;

      if(likePlot == TRUE){
        graphics::par();
        stats::ts.plot(mcmc$ll.trace, xlab = "Iteration", ylab = "log-likelihood",
                       main = paste("subset ", i, sep=""));
      }
    }
  }

  cat(paste("Time elapsed: ", round(as.numeric(proc.time()[1] - ptime) / 60, 2),
            " minutes", sep = ""), "\n");

  out[["info"]] = list(p = p, l = l);

  class(out) = "psds"; # Assign S3 class to object

  return(out);

}
