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
