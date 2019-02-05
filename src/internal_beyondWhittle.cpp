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
