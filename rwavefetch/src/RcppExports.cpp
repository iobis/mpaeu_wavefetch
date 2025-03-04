// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// meanC
double meanC(NumericVector x);
RcppExport SEXP _rwavefetch_meanC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(meanC(x));
    return rcpp_result_gen;
END_RCPP
}
// isitcoast
IntegerMatrix isitcoast(IntegerMatrix land);
RcppExport SEXP _rwavefetch_isitcoast(SEXP landSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type land(landSEXP);
    rcpp_result_gen = Rcpp::wrap(isitcoast(land));
    return rcpp_result_gen;
END_RCPP
}
// isitnearcoast
IntegerMatrix isitnearcoast(IntegerMatrix land, int dx);
RcppExport SEXP _rwavefetch_isitnearcoast(SEXP landSEXP, SEXP dxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type land(landSEXP);
    Rcpp::traits::input_parameter< int >::type dx(dxSEXP);
    rcpp_result_gen = Rcpp::wrap(isitnearcoast(land, dx));
    return rcpp_result_gen;
END_RCPP
}
// coastal_wave_fetch
NumericMatrix coastal_wave_fetch(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit, double land1llx, double land1uly, double land1cellszx, double land1cellszy, double land2llx, double land2uly, double land2cellszx, double land2cellszy, int verbose);
RcppExport SEXP _rwavefetch_coastal_wave_fetch(SEXP land1SEXP, SEXP land2SEXP, SEXP dxSEXP, SEXP dwxSEXP, SEXP clSEXP, SEXP jitSEXP, SEXP land1llxSEXP, SEXP land1ulySEXP, SEXP land1cellszxSEXP, SEXP land1cellszySEXP, SEXP land2llxSEXP, SEXP land2ulySEXP, SEXP land2cellszxSEXP, SEXP land2cellszySEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type land1(land1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type land2(land2SEXP);
    Rcpp::traits::input_parameter< int >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< int >::type dwx(dwxSEXP);
    Rcpp::traits::input_parameter< int >::type cl(clSEXP);
    Rcpp::traits::input_parameter< int >::type jit(jitSEXP);
    Rcpp::traits::input_parameter< double >::type land1llx(land1llxSEXP);
    Rcpp::traits::input_parameter< double >::type land1uly(land1ulySEXP);
    Rcpp::traits::input_parameter< double >::type land1cellszx(land1cellszxSEXP);
    Rcpp::traits::input_parameter< double >::type land1cellszy(land1cellszySEXP);
    Rcpp::traits::input_parameter< double >::type land2llx(land2llxSEXP);
    Rcpp::traits::input_parameter< double >::type land2uly(land2ulySEXP);
    Rcpp::traits::input_parameter< double >::type land2cellszx(land2cellszxSEXP);
    Rcpp::traits::input_parameter< double >::type land2cellszy(land2cellszySEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(coastal_wave_fetch(land1, land2, dx, dwx, cl, jit, land1llx, land1uly, land1cellszx, land1cellszy, land2llx, land2uly, land2cellszx, land2cellszy, verbose));
    return rcpp_result_gen;
END_RCPP
}
// coastal_wave_orientation
NumericMatrix coastal_wave_orientation(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit, double land1llx, double land1uly, double land1cellszx, double land1cellszy, double land2llx, double land2uly, double land2cellszx, double land2cellszy, int verbose);
RcppExport SEXP _rwavefetch_coastal_wave_orientation(SEXP land1SEXP, SEXP land2SEXP, SEXP dxSEXP, SEXP dwxSEXP, SEXP clSEXP, SEXP jitSEXP, SEXP land1llxSEXP, SEXP land1ulySEXP, SEXP land1cellszxSEXP, SEXP land1cellszySEXP, SEXP land2llxSEXP, SEXP land2ulySEXP, SEXP land2cellszxSEXP, SEXP land2cellszySEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type land1(land1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type land2(land2SEXP);
    Rcpp::traits::input_parameter< int >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< int >::type dwx(dwxSEXP);
    Rcpp::traits::input_parameter< int >::type cl(clSEXP);
    Rcpp::traits::input_parameter< int >::type jit(jitSEXP);
    Rcpp::traits::input_parameter< double >::type land1llx(land1llxSEXP);
    Rcpp::traits::input_parameter< double >::type land1uly(land1ulySEXP);
    Rcpp::traits::input_parameter< double >::type land1cellszx(land1cellszxSEXP);
    Rcpp::traits::input_parameter< double >::type land1cellszy(land1cellszySEXP);
    Rcpp::traits::input_parameter< double >::type land2llx(land2llxSEXP);
    Rcpp::traits::input_parameter< double >::type land2uly(land2ulySEXP);
    Rcpp::traits::input_parameter< double >::type land2cellszx(land2cellszxSEXP);
    Rcpp::traits::input_parameter< double >::type land2cellszy(land2cellszySEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(coastal_wave_orientation(land1, land2, dx, dwx, cl, jit, land1llx, land1uly, land1cellszx, land1cellszy, land2llx, land2uly, land2cellszx, land2cellszy, verbose));
    return rcpp_result_gen;
END_RCPP
}
// coastal_wave_direction
NumericMatrix coastal_wave_direction(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit, double land1llx, double land1uly, double land1cellszx, double land1cellszy, double land2llx, double land2uly, double land2cellszx, double land2cellszy, int verbose);
RcppExport SEXP _rwavefetch_coastal_wave_direction(SEXP land1SEXP, SEXP land2SEXP, SEXP dxSEXP, SEXP dwxSEXP, SEXP clSEXP, SEXP jitSEXP, SEXP land1llxSEXP, SEXP land1ulySEXP, SEXP land1cellszxSEXP, SEXP land1cellszySEXP, SEXP land2llxSEXP, SEXP land2ulySEXP, SEXP land2cellszxSEXP, SEXP land2cellszySEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type land1(land1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type land2(land2SEXP);
    Rcpp::traits::input_parameter< int >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< int >::type dwx(dwxSEXP);
    Rcpp::traits::input_parameter< int >::type cl(clSEXP);
    Rcpp::traits::input_parameter< int >::type jit(jitSEXP);
    Rcpp::traits::input_parameter< double >::type land1llx(land1llxSEXP);
    Rcpp::traits::input_parameter< double >::type land1uly(land1ulySEXP);
    Rcpp::traits::input_parameter< double >::type land1cellszx(land1cellszxSEXP);
    Rcpp::traits::input_parameter< double >::type land1cellszy(land1cellszySEXP);
    Rcpp::traits::input_parameter< double >::type land2llx(land2llxSEXP);
    Rcpp::traits::input_parameter< double >::type land2uly(land2ulySEXP);
    Rcpp::traits::input_parameter< double >::type land2cellszx(land2cellszxSEXP);
    Rcpp::traits::input_parameter< double >::type land2cellszy(land2cellszySEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(coastal_wave_direction(land1, land2, dx, dwx, cl, jit, land1llx, land1uly, land1cellszx, land1cellszy, land2llx, land2uly, land2cellszx, land2cellszy, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rwavefetch_meanC", (DL_FUNC) &_rwavefetch_meanC, 1},
    {"_rwavefetch_isitcoast", (DL_FUNC) &_rwavefetch_isitcoast, 1},
    {"_rwavefetch_isitnearcoast", (DL_FUNC) &_rwavefetch_isitnearcoast, 2},
    {"_rwavefetch_coastal_wave_fetch", (DL_FUNC) &_rwavefetch_coastal_wave_fetch, 15},
    {"_rwavefetch_coastal_wave_orientation", (DL_FUNC) &_rwavefetch_coastal_wave_orientation, 15},
    {"_rwavefetch_coastal_wave_direction", (DL_FUNC) &_rwavefetch_coastal_wave_direction, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_rwavefetch(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
