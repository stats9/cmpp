// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Initialize
void Initialize(NumericMatrix features, NumericVector x, IntegerVector delta1, IntegerVector delta2, double h);
RcppExport SEXP _cmpp_Initialize(SEXP featuresSEXP, SEXP xSEXP, SEXP delta1SEXP, SEXP delta2SEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type features(featuresSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type delta2(delta2SEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Initialize(features, x, delta1, delta2, h);
    return R_NilValue;
END_RCPP
}
// cdf_gomp
double cdf_gomp(double x, double alpha, double beta);
RcppExport SEXP _cmpp_cdf_gomp(SEXP xSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(cdf_gomp(x, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// pdf_gomp
double pdf_gomp(double x, double alpha, double beta);
RcppExport SEXP _cmpp_pdf_gomp(SEXP xSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(pdf_gomp(x, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// GetDim
Rcpp::List GetDim();
RcppExport SEXP _cmpp_GetDim() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(GetDim());
    return rcpp_result_gen;
END_RCPP
}
// LogLike1
SEXP LogLike1(SEXP param);
RcppExport SEXP _cmpp_LogLike1(SEXP paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type param(paramSEXP);
    rcpp_result_gen = Rcpp::wrap(LogLike1(param));
    return rcpp_result_gen;
END_RCPP
}
// compute_grad
SEXP compute_grad(SEXP param);
RcppExport SEXP _cmpp_compute_grad(SEXP paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type param(paramSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_grad(param));
    return rcpp_result_gen;
END_RCPP
}
// compute_hessian
SEXP compute_hessian(SEXP param);
RcppExport SEXP _cmpp_compute_hessian(SEXP paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type param(paramSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_hessian(param));
    return rcpp_result_gen;
END_RCPP
}
// makeMat
Eigen::MatrixXd makeMat(int n, int m, double value);
RcppExport SEXP _cmpp_makeMat(SEXP nSEXP, SEXP mSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(makeMat(n, m, value));
    return rcpp_result_gen;
END_RCPP
}
// Cleanup
void Cleanup();
RcppExport SEXP _cmpp_Cleanup() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Cleanup();
    return R_NilValue;
END_RCPP
}
// bootstrap_variance
List bootstrap_variance(NumericMatrix features, NumericVector x, IntegerVector delta1, IntegerVector delta2, NumericVector initial_params, int n_bootstrap, std::string optimMethod);
RcppExport SEXP _cmpp_bootstrap_variance(SEXP featuresSEXP, SEXP xSEXP, SEXP delta1SEXP, SEXP delta2SEXP, SEXP initial_paramsSEXP, SEXP n_bootstrapSEXP, SEXP optimMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type features(featuresSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type delta2(delta2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initial_params(initial_paramsSEXP);
    Rcpp::traits::input_parameter< int >::type n_bootstrap(n_bootstrapSEXP);
    Rcpp::traits::input_parameter< std::string >::type optimMethod(optimMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap_variance(features, x, delta1, delta2, initial_params, n_bootstrap, optimMethod));
    return rcpp_result_gen;
END_RCPP
}
// F_cdf_rcpp
double F_cdf_rcpp(NumericVector Params, NumericVector Z, double x);
RcppExport SEXP _cmpp_F_cdf_rcpp(SEXP ParamsSEXP, SEXP ZSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(F_cdf_rcpp(Params, Z, x));
    return rcpp_result_gen;
END_RCPP
}
// f_pdf_rcpp
double f_pdf_rcpp(NumericVector Params, NumericVector Z, double x);
RcppExport SEXP _cmpp_f_pdf_rcpp(SEXP ParamsSEXP, SEXP ZSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(f_pdf_rcpp(Params, Z, x));
    return rcpp_result_gen;
END_RCPP
}
// log_f_rcpp
double log_f_rcpp(NumericVector Params);
RcppExport SEXP _cmpp_log_f_rcpp(SEXP ParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(log_f_rcpp(Params));
    return rcpp_result_gen;
END_RCPP
}
// compute_log_f_gradient_rcpp
NumericVector compute_log_f_gradient_rcpp(NumericVector Params);
RcppExport SEXP _cmpp_compute_log_f_gradient_rcpp(SEXP ParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log_f_gradient_rcpp(Params));
    return rcpp_result_gen;
END_RCPP
}
// compute_log_f_hessian_rcpp
NumericMatrix compute_log_f_hessian_rcpp(NumericVector Params);
RcppExport SEXP _cmpp_compute_log_f_hessian_rcpp(SEXP ParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log_f_hessian_rcpp(Params));
    return rcpp_result_gen;
END_RCPP
}
// F_cdf_rcpp2
double F_cdf_rcpp2(NumericVector Params, NumericVector Z, double x);
RcppExport SEXP _cmpp_F_cdf_rcpp2(SEXP ParamsSEXP, SEXP ZSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(F_cdf_rcpp2(Params, Z, x));
    return rcpp_result_gen;
END_RCPP
}
// f_pdf_rcpp2
double f_pdf_rcpp2(NumericVector Params, NumericVector Z, double x);
RcppExport SEXP _cmpp_f_pdf_rcpp2(SEXP ParamsSEXP, SEXP ZSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(f_pdf_rcpp2(Params, Z, x));
    return rcpp_result_gen;
END_RCPP
}
// log_f_rcpp2
double log_f_rcpp2(NumericVector Params);
RcppExport SEXP _cmpp_log_f_rcpp2(SEXP ParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(log_f_rcpp2(Params));
    return rcpp_result_gen;
END_RCPP
}
// compute_log_f_gradient_rcpp2
NumericVector compute_log_f_gradient_rcpp2(NumericVector Params);
RcppExport SEXP _cmpp_compute_log_f_gradient_rcpp2(SEXP ParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log_f_gradient_rcpp2(Params));
    return rcpp_result_gen;
END_RCPP
}
// F_cdf_rcpp3
double F_cdf_rcpp3(NumericVector Params, NumericVector Z, double x);
RcppExport SEXP _cmpp_F_cdf_rcpp3(SEXP ParamsSEXP, SEXP ZSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(F_cdf_rcpp3(Params, Z, x));
    return rcpp_result_gen;
END_RCPP
}
// f_pdf_rcpp3
double f_pdf_rcpp3(NumericVector Params, NumericVector Z, double x);
RcppExport SEXP _cmpp_f_pdf_rcpp3(SEXP ParamsSEXP, SEXP ZSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(f_pdf_rcpp3(Params, Z, x));
    return rcpp_result_gen;
END_RCPP
}
// log_f_rcpp3
double log_f_rcpp3(NumericVector Params);
RcppExport SEXP _cmpp_log_f_rcpp3(SEXP ParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(log_f_rcpp3(Params));
    return rcpp_result_gen;
END_RCPP
}
// compute_log_f_gradient_rcpp3
NumericVector compute_log_f_gradient_rcpp3(NumericVector Params);
RcppExport SEXP _cmpp_compute_log_f_gradient_rcpp3(SEXP ParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Params(ParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log_f_gradient_rcpp3(Params));
    return rcpp_result_gen;
END_RCPP
}
// GetData
Rcpp::List GetData();
RcppExport SEXP _cmpp_GetData() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(GetData());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cmpp_Initialize", (DL_FUNC) &_cmpp_Initialize, 5},
    {"_cmpp_cdf_gomp", (DL_FUNC) &_cmpp_cdf_gomp, 3},
    {"_cmpp_pdf_gomp", (DL_FUNC) &_cmpp_pdf_gomp, 3},
    {"_cmpp_GetDim", (DL_FUNC) &_cmpp_GetDim, 0},
    {"_cmpp_LogLike1", (DL_FUNC) &_cmpp_LogLike1, 1},
    {"_cmpp_compute_grad", (DL_FUNC) &_cmpp_compute_grad, 1},
    {"_cmpp_compute_hessian", (DL_FUNC) &_cmpp_compute_hessian, 1},
    {"_cmpp_makeMat", (DL_FUNC) &_cmpp_makeMat, 3},
    {"_cmpp_Cleanup", (DL_FUNC) &_cmpp_Cleanup, 0},
    {"_cmpp_bootstrap_variance", (DL_FUNC) &_cmpp_bootstrap_variance, 7},
    {"_cmpp_F_cdf_rcpp", (DL_FUNC) &_cmpp_F_cdf_rcpp, 3},
    {"_cmpp_f_pdf_rcpp", (DL_FUNC) &_cmpp_f_pdf_rcpp, 3},
    {"_cmpp_log_f_rcpp", (DL_FUNC) &_cmpp_log_f_rcpp, 1},
    {"_cmpp_compute_log_f_gradient_rcpp", (DL_FUNC) &_cmpp_compute_log_f_gradient_rcpp, 1},
    {"_cmpp_compute_log_f_hessian_rcpp", (DL_FUNC) &_cmpp_compute_log_f_hessian_rcpp, 1},
    {"_cmpp_F_cdf_rcpp2", (DL_FUNC) &_cmpp_F_cdf_rcpp2, 3},
    {"_cmpp_f_pdf_rcpp2", (DL_FUNC) &_cmpp_f_pdf_rcpp2, 3},
    {"_cmpp_log_f_rcpp2", (DL_FUNC) &_cmpp_log_f_rcpp2, 1},
    {"_cmpp_compute_log_f_gradient_rcpp2", (DL_FUNC) &_cmpp_compute_log_f_gradient_rcpp2, 1},
    {"_cmpp_F_cdf_rcpp3", (DL_FUNC) &_cmpp_F_cdf_rcpp3, 3},
    {"_cmpp_f_pdf_rcpp3", (DL_FUNC) &_cmpp_f_pdf_rcpp3, 3},
    {"_cmpp_log_f_rcpp3", (DL_FUNC) &_cmpp_log_f_rcpp3, 1},
    {"_cmpp_compute_log_f_gradient_rcpp3", (DL_FUNC) &_cmpp_compute_log_f_gradient_rcpp3, 1},
    {"_cmpp_GetData", (DL_FUNC) &_cmpp_GetData, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_cmpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
