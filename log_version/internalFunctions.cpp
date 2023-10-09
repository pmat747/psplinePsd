#include <Rcpp.h>
using namespace Rcpp;

//' C++ function for calculating the sum of the log of two log-values,i.e., logplus(log(a),log(b))=log(a+b).
//' @keywords internal
// [[Rcpp::export]]
double logplus(double x, double y){

  if(x>y){

    return x + log(1.0+exp(y-x));

  }else{

    return y + log(1.0+exp(x-y));

  }
}

//' C++ function for ...
//' @keywords internal
// [[Rcpp::export]]
double logplusvec( NumericVector & x){

  int n = x.size();

  double r = -DBL_MAX;

  for(int i = 0; i < n; i++){

    r = logplus(r, x[i]);

  }

  return r;

}
