#ifndef MYFUNCTIONS_H
#define MYFUNCTIONS_H

#include <Rcpp.h>

extern "C" SEXP cdf_gomp(SEXP x, SEXP alpha, SEXP beta);
extern "C" SEXP compute_grad(SEXP Param);
extern "C" SEXP GrLike(SEXP Param, SEXP grad);
extern "C" SEXP Initialize(SEXP features, SEXP x, SEXP delta1, SEXP delta2, SEXP h);
extern "C" SEXP LogLike1(SEXP Param);
extern "C" SEXP pdf_gomp(SEXP x, SEXP alpha, SEXP beta);

#endif // MYFUNCTIONS_H
