#include <Rcpp.h>
using namespace Rcpp;

//' C++ function for building a density mixture, given mixture weights and functions.
//' @keywords internal
// [[Rcpp::export]]
NumericVector densityMixture(NumericVector weights, NumericMatrix densities) {
  if (weights.size() != densities.nrow()) {
    return(NumericVector());
  }
  const unsigned n = densities.ncol();
  NumericVector res(n);
  for (unsigned omega = 0; omega < n; ++omega) {
    res[omega] = 0.0;
  }
  for (unsigned j = 0; j < weights.size(); ++j) {
    for (unsigned omega = 0; omega < n; ++omega) {
      res[omega] += weights[j] * densities(j, omega);
    }
  }
  return(res);
}

//' C++ help function to redundantly roll out a PSD to length n
//' @keywords internal
// [[Rcpp::export]]
NumericVector unrollPsd(NumericVector qPsd, unsigned n) {
  NumericVector q(n);
  q[0] = qPsd[0];
  const unsigned N = (n-1)/2;
  for (unsigned i = 1; i <= N; ++i) {
    const unsigned j = 2 * i - 1;
    q[j] = qPsd[i];
    q[j+1] = qPsd[i];
  }
  if (!(n % 2)) {
    q[n-1] = qPsd[qPsd.size() - 1];
  }
  return(q);
}
