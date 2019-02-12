#' Difference Matrix. It works for second & third order penalties
#' @keywords internal
diffMatrix = function(k, d = 2){

  if(k<1)stop("k must be greater than 0")

  if(d == 1){

    x = c(-1, 1, rep(c(rep(0, k-1), -1, 1), k-2));
    x = matrix(x, ncol = k, byrow = TRUE);
    return(x)

  }else if(d == 2){

    x = c(1,-2,1, rep(c(rep(0, k-2), 1,-2,1), k-3));
    x = matrix(x, ncol = k, byrow = TRUE);
    return(x);

  }else if(d == 3){

    x = c(1, -3, 3, -1, rep(c(rep(0, k-3), 1,-3,3,-1), k-4));
    x = matrix(x, ncol = k, byrow = TRUE);
    return(x);

  }else{

    stop("diffMatrix has been defined for d = 1, 2 and 3");

  }
}
