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
knotLoc = function(data, k, degree, eqSpaced = FALSE, nfreqbin = 1, 
                   wfreqbin = NULL){
  
  # data     = numeric vector containing the data
  # k        = number of B-spline densities
  # degree   = positive integer specifying the degree of the B-spline densities
  # eqSpaced = If it is TRUE, it generates equidistant knots.
  # nfreqbin = number of partitions in the frequency domain or
  #            a vector containing the location of the partitions between ]0,pi[.
  # wfreqbin = percentage of knots assigned to each bin.  If it is not specified
  #            equal probabilities are assigned.
  
  K = k - degree + 1; # number of internal knots [0,1]
  
  if(eqSpaced == TRUE){
    
    knots = seq(from = 0, to = 1, length = K); # Equidistant knots
    
    return(knots);
    
  }
  
  # Periodogram
  data  = data - mean(data);
  FZ    = fast_ft(data); # FFT data to frequency domain.  NOTE: Must be mean-centred.
  pdgrm = abs(FZ) ^ 2; # Periodogram: NOTE: the length is n here.
  #pdgrm = log(pdgrm);  # In case using log(pdgrm)
  #pdgrm[1] = pdgrm[2]; # In case using log(pdgrm)
  N     = length(pdgrm);
  
  # Defining and checking "wfreqbin"
  if(length(nfreqbin) == 1){ 
    
    # 'nfreqbin' defined as the number of frequency bins
    
    if(is.null(wfreqbin)){
      
      wfreqbin = rep(1, nfreqbin) / nfreqbin; # Equal weight for each bin
      
    }else{
      
      wfreqbin = wfreqbin / sum(wfreqbin); # normalising the weights
      
    }
    
    # Positions in pdgrm to be considered in each frequency bin
    index = round(seq(1, N, length = nfreqbin+1), 0);
    index = cbind(c(0, index[-c(1, length(index))])+1, index[-1]);
    
    Nindex = ifelse(is.null(dim(index)[1]), 1, dim(index)[1]);
    
    eqval = seq(0, 1, length = N);
    eqval = cbind(eqval[index[,1]], eqval[index[,2]]);# Each row represents the bin
    
  }else{
    
    # 'nfreqbin' defined as vector
    
    if(is.null(wfreqbin)){
      
      n_wfreqbin = length(nfreqbin) + 1; # number of elements in wfreqbin
      
      wfreqbin = rep(1, n_wfreqbin) / n_wfreqbin;# Equal weight for each bin
      
    }else{
      
      if(length(nfreqbin) + 1 != length(wfreqbin)){
        stop("The length of 'wfreqbin' must be equal to the length of 'nfreqbin + 1'");
      }
      
      wfreqbin = wfreqbin / sum(wfreqbin); # normalising the weights
      n_wfreqbin = length(wfreqbin);
      
    }
    
    nfreqbin = nfreqbin[order(nfreqbin)]; # Ordering (increasingly)
    eqval = c(0, nfreqbin / pi, 1); # Transforming to interval [0,1]
    eqval = cbind(eqval[-length(eqval)], eqval[-1]); # Each row represents the bin
    
    j = seq(0, 1, length = N);
    s = seq(1,N);
    index = NULL;
    
    for(i in 1:n_wfreqbin){
      
      cond = (j >= eqval[i, 1]) & (j <= eqval[i,2]);
      
      index = rbind(index, c(min(s[cond]), max(s[cond])));
    }
    
    Nindex = dim(index)[1];# Each row represents the bin position in pdgrm
    
  }
  
  K    = K - 2; # to include 0 and 1 in the knot vector
  kvec = round(K * wfreqbin, 0);

  knots = NULL;
  
  for(i in 1:Nindex){
    
    aux_pdgrm = pdgrm[index[i,1]:index[i,2]];
    
    aux  = sqrt(aux_pdgrm);
    #aux  = aux_pdgrm; # In case using log(pdgrm)
    dens = abs(aux - mean(aux))/stats::sd(aux);
    
    Naux = length(aux);
    
    dens = dens / sum(dens);
    cumf = cumsum(dens);
    
    # Distribution function 
    df  = stats::approxfun(x = seq(eqval[i,1],eqval[i,2], length = Naux),
                           y = cumf, yleft = 0, yright = 1);
    
    # Inverse distribution 
    vec   = seq(eqval[i,1], eqval[i,2], length = Naux);
    dfvec = df(vec);
    
    invDf = stats::approxfun(x = dfvec, y = vec, yleft = NA, yright = NA);
    
    v     = seq(0, 1, length = kvec[i] + 2);
    v     = v[-c(1, kvec[i] + 2)];
    
    knots = c(knots, invDf(v));
    
  }
  
  knots = c(0, knots, 1);
  
  return(knots);
  
}
