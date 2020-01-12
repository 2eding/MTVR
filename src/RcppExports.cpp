// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// generateCovBand
void generateCovBand(int windowSize, std::string corrMatrix, std::string output);
RcppExport SEXP _MTVR_generateCovBand(SEXP windowSizeSEXP, SEXP corrMatrixSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type windowSize(windowSizeSEXP);
    Rcpp::traits::input_parameter< std::string >::type corrMatrix(corrMatrixSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    generateCovBand(windowSize, corrMatrix, output);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MTVR_generateCovBand", (DL_FUNC) &_MTVR_generateCovBand, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_MTVR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
