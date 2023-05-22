// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// logplus
double logplus(double x, double y);
RcppExport SEXP _psplinePsd_logplus(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(logplus(x, y));
    return rcpp_result_gen;
END_RCPP
}
// logplusvec
double logplusvec(NumericVector& x);
RcppExport SEXP _psplinePsd_logplusvec(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logplusvec(x));
    return rcpp_result_gen;
END_RCPP
}
// densityMixture
NumericVector densityMixture(NumericVector weights, NumericMatrix densities);
RcppExport SEXP _psplinePsd_densityMixture(SEXP weightsSEXP, SEXP densitiesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type densities(densitiesSEXP);
    rcpp_result_gen = Rcpp::wrap(densityMixture(weights, densities));
    return rcpp_result_gen;
END_RCPP
}
// unrollPsd
NumericVector unrollPsd(NumericVector qPsd, unsigned n);
RcppExport SEXP _psplinePsd_unrollPsd(SEXP qPsdSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type qPsd(qPsdSEXP);
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(unrollPsd(qPsd, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_psplinePsd_logplus", (DL_FUNC) &_psplinePsd_logplus, 2},
    {"_psplinePsd_logplusvec", (DL_FUNC) &_psplinePsd_logplusvec, 1},
    {"_psplinePsd_densityMixture", (DL_FUNC) &_psplinePsd_densityMixture, 2},
    {"_psplinePsd_unrollPsd", (DL_FUNC) &_psplinePsd_unrollPsd, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_psplinePsd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
