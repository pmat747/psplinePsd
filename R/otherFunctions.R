#' diffMatrix calculates the difference penalty matrix
#' The third order difference penalty matrix is difference from the one
#' produced by the other functions.  However, the penalty matrix t(D)%*%D
#' remains the same
#' @keywords internal
diffMatrix = function(k, d = 2){

  if( (d<1) || (d %% 1 != 0) )stop("d must be a positive integer value");
  if( (k<1) || (k %% 1 != 0) )stop("k must be a positive integer value");
  if(d >= k)stop("d must be lower than k");

  out = diag(k);

  for(i in 1:d){

    out = diff(out);

  }
  return(out)
}

#' This function produces a matrix which each row contains the first and last indexes
#'  to split a time series of length "n"
#'
#' n  = time series length
#' l  = interval length
#' p  = overlapping percentage
#' eq = last interval has length "l"
#' @keywords internal
ints = function(n, l, p = 00, eq = TRUE){

  if( (n<=0) || (n %% 1 != 0)) stop("n must be a positive integer")
  if( (l<=0) || (l %% 1 != 0)) stop("l must be a positive integer")
  if( (p<0) || (p>=100)) stop("p must be a positive number between 0 and 100")
  if( l >= n ) stop("l must be lower than n")

  # This version yields eve interval lengths
  if (l %% 2 != 0) stop("this version of bsplinePsd must have l even")

  ovp  = round(l*p/100,0); # number of overlaped points
  a    = l - ovp + 1;
  col1 = c(1, seq(from = a, to = n - l + a -1, by = a-1));

  index = cbind(col1, col1 + l-1);
  # checking the lengths
  # index[,2] - index[,1] + 1 == l
  # diff(index[,1]) == a-1
  # diff(index[,2]) == a-1

  # fixing last interval
  index[dim(index)[1], 2] = n;

  colnames(index) = NULL;
  rownames(index) = NULL;

  cat(paste("The number of data subsets is ",dim(index)[1], sep=""),"\n");

  if(eq == TRUE){

    lint = dim(index)[1]
    index[lint, 1] = n - l + 1;
    # index[,2] - index[,1] == l

    x = index[lint-1, 2] - index[lint, 1] + 1;
    x = round(100 * x / l, 2);

    if(x!=p){
      cat(paste("Last interval overlaps ", x, "%", sep = ""), "\n");
    }

  }else{

    aux = index[dim(index)[1], 2] - index[dim(index)[1], 1] + 1;

    if(aux %% 2 != 0){

      # last interval is been made even
      index[dim(index)[1], 1] = index[dim(index)[1], 1] - 1;
      aux = index[dim(index)[1], 2] - index[dim(index)[1], 1] + 1;
    }

    if(aux != l){
      cat(paste("Last interval contains ", aux, " observations", sep = ""), "\n");
    }
  }

  return(index);
}


#' This function produces the knot location
#'
#' data     = data
#' k        = number of B-spline densities
#' degree   = positive integer specifying the degree of the B-spline densities (default is 3)
#' eqSpaced = if TRUE, it produces equally spaced knots, otherwise it allocates the knots according to the periodogram
#' @importFrom stats approxfun sd
#' @keywords internal
knotLoc = function(data, k, degree, eqSpaced = FALSE){

  if(eqSpaced == TRUE){

    knots = seq(from = 0, to = 1, length = k - degree + 1);

    return(knots);

  }

  data  = data - mean(data);
  FZ    = fast_ft(data); # FFT data to frequency domain.  NOTE: Must be mean-centred.
  pdgrm = abs(FZ) ^ 2; # Periodogram: NOTE: the length is n here.

  #aux   = (pdgrm - mean(pdgrm)) / sd(pdgrm);
  aux   = sqrt(pdgrm)
  dens  = abs(aux - mean(aux))/stats::sd(aux)

  N     = length(pdgrm);

  # function based on peridogram
  #f = stats::approxfun(x = seq(0,1,length = N), y = abs(aux)/sum(abs(aux)),
  #c                     yleft = 0, yright = 0);

  # cumulative values of f function
  #n    = 1000;
  #dens = f(seq(from = 0, to = 1, length = n));
  dens = dens / sum(dens);
  cumf = cumsum(dens);

  # distribution function of f
  df  = stats::approxfun(x = seq(0,1,length = N),
                         y = cumf, yleft = 0, yright = 1);

  # inverse distribution of f
  invDf = stats::approxfun(x = df(seq(0,1,length = N)),
                           y = seq(0,1,length = N), yleft = 0, yright = 1);

  # equally spaced knots
  knots = seq(0,1,length = k - degree + 1);

  # knots based on peaks of the periodogram
  knots = invDf(knots);

  return(knots);

}
