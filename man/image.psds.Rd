% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_psdsPlot_s3.R
\name{image.psds}
\alias{image.psds}
\title{Plot method for psds class}
\usage{
\method{image}{psds}(x, func = "median", logZ = TRUE, zoomFreq = c(0, 1), fs = 16384, ...)
}
\arguments{
\item{x}{an object of class psds}

\item{func}{statistics to be applied to \code{x}, e.g., \code{median} or \code{mean}}

\item{logZ}{logical value (default is TRUE) to determine if the PSD should be on log scale}

\item{zoomFreq}{interval which specifies the frequency axis to be plotted, for instance, c(0.2, 0.7) leaves out the 20\% and 30\% of the low and high frequencies, respectively.  The default value is c(0,1), i.e., all the frequencies}

\item{fs}{sampling frequency}

\item{...}{other graphical parameters from the image.plot function}
}
\value{
plot of the estimate of the (log) spectrum
}
\description{
This function produces a spectogram from a psds object.
}
\examples{
\dontrun{

set.seed(1)

# Generate AR(1) data with rho = 0.9
n = 128
data = arima.sim(n, model = list(ar = 0.9));
data = data - mean(data);

# Spectrum estimate via p-splines
spec = spec_pspline(data, l=30, p=20, Ntotal1=3000, burnin1=1000, thin1=10,
                    Ntotal=3000, burnin=1000, thin=10, k=30);

image(spec) # Plot log PSD (see documentation of iplot.psds)
}
}
\seealso{
\link{gibbs_pspline}
}
