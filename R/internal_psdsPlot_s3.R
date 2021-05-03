#' Plot method for psds class
#' @description This function produces a spectogram from a psds object.
#' @export
#' @param x an object of class psds
#' @param func statistics to be applied to \code{x}, e.g., \code{median} or \code{mean}
#' @param logZ logical value (default is TRUE) to determine if the PSD should be on log scale
#' @param zoomFreq interval which specifies the frequency axis to be plotted, for instance, c(0.2, 0.7) leaves out the 20\% and 30\% of the low and high frequencies, respectively.  The default value is c(0,1), i.e., all the frequencies
#' @param fs sampling frequency
#' @param ... other graphical parameters from the image.plot function
#' @return plot of the estimate of the (log) spectrum
#' @seealso \link{gibbs_pspline}
#' @importFrom fields image.plot
#'
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
#' # Spectrum estimate via p-splines
#' spec = spec_pspline(data, l=30, p=20, Ntotal1=3000, burnin1=1000, thin1=10,
#'                     Ntotal=3000, burnin=1000, thin=10, k=30);
#'
#' image(spec) # Plot log PSD (see documentation of iplot.psds)
#' }
image.psds = function(x, func = 'median', logZ = TRUE, zoomFreq = c(0,1), fs = 16384, ...) {  # Plot method for "psds" class

  if((zoomFreq[1]<0) || (zoomFreq[2]>1) || (zoomFreq[1] > zoomFreq[2])){
    stop("zoomFreq must be a vector c(a,b) with values 0 <= a < b <= 1")
  }

  N   = length(x$psds); # number of fpsd.sample elements
  po  = x$info$p/100;
  pc  = 1-po;
  aux = NULL;

  for(i in 1:N){
    aux = cbind(aux, apply(x$psds[[i]], 1, get(func))); # MEDIAN
  }

  aux = sqrt(aux);

  if(logZ == TRUE){

    aux = log10(aux);

  }

  i <- 1:N;
  x <- round(x$info$l / fs * (i-1) * pc * 1000, 2); # number of intervals
  #y <- seq(0,1,length=dim(aux)[1]) # number of frequencies
  y <- seq(1,fs/2,length=dim(aux)[1]);

  if(all(zoomFreq == c(0,1))){

    z <- t(aux);
    z <- z[, -c(1, dim(z)[2])]; #deleting first and last row in image
    y <- y[-c(1, dim(z)[2])]; # adjustment

  }else{

    index    = round(dim(aux)[1] * zoomFreq);
    index    = seq(index[1] + 1, index[2]);

    y <- y[index]
    z <- t(aux[index, ])

  }

  fields::image.plot(x, y, z, xlab = "Time [ms]", ylab = "Frequency [Hz]",
                     main = paste(func), ...); # yaxt='n'

}
