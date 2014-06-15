#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#if !defined(__MATHUTILS_H)
#define __MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

double rsCorrelation(const double* X, const double *Y, const size_t length);
double rsZCorrelation(const double* X, const double* Y, const size_t length);
double rsFastZCorrelation(const double* X, const double* Y, const size_t length);
double rsTCorrelation(const double* X, const double* Y, const size_t length);
double rsMonteCarloZCorrelation(const double* X, const double* Y, const size_t length, const unsigned int repetitions, const unsigned int samplingSize);

#ifdef __cplusplus
}
#endif

#endif
