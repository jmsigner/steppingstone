// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// tpm_func
NumericMatrix tpm_func(double alpha, NumericVector omegas, NumericMatrix resources, NumericVector d, int nc);
RcppExport SEXP steppingstone_tpm_func(SEXP alphaSEXP, SEXP omegasSEXP, SEXP resourcesSEXP, SEXP dSEXP, SEXP ncSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type omegas(omegasSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type resources(resourcesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    rcpp_result_gen = Rcpp::wrap(tpm_func(alpha, omegas, resources, d, nc));
    return rcpp_result_gen;
END_RCPP
}
// walk
List walk(NumericMatrix tpm, int n, NumericVector xy0, int nc, NumericVector dp, int init_dir, int boundary, int max_try);
RcppExport SEXP steppingstone_walk(SEXP tpmSEXP, SEXP nSEXP, SEXP xy0SEXP, SEXP ncSEXP, SEXP dpSEXP, SEXP init_dirSEXP, SEXP boundarySEXP, SEXP max_trySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type tpm(tpmSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xy0(xy0SEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dp(dpSEXP);
    Rcpp::traits::input_parameter< int >::type init_dir(init_dirSEXP);
    Rcpp::traits::input_parameter< int >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< int >::type max_try(max_trySEXP);
    rcpp_result_gen = Rcpp::wrap(walk(tpm, n, xy0, nc, dp, init_dir, boundary, max_try));
    return rcpp_result_gen;
END_RCPP
}
// ud_func
List ud_func(NumericMatrix tpm, int n, NumericVector xy0, int nc, int burnin, NumericVector dp, int boundary, int init_dir, int max_try);
RcppExport SEXP steppingstone_ud_func(SEXP tpmSEXP, SEXP nSEXP, SEXP xy0SEXP, SEXP ncSEXP, SEXP burninSEXP, SEXP dpSEXP, SEXP boundarySEXP, SEXP init_dirSEXP, SEXP max_trySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type tpm(tpmSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xy0(xy0SEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dp(dpSEXP);
    Rcpp::traits::input_parameter< int >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< int >::type init_dir(init_dirSEXP);
    Rcpp::traits::input_parameter< int >::type max_try(max_trySEXP);
    rcpp_result_gen = Rcpp::wrap(ud_func(tpm, n, xy0, nc, burnin, dp, boundary, init_dir, max_try));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"steppingstone_tpm_func", (DL_FUNC) &steppingstone_tpm_func, 5},
    {"steppingstone_walk", (DL_FUNC) &steppingstone_walk, 8},
    {"steppingstone_ud_func", (DL_FUNC) &steppingstone_ud_func, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_steppingstone(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
