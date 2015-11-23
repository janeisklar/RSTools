#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>

#if !defined(__STATS_MATHUTILS_H)
#define __STATS_MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif
    
float  rsFastErfInv(float x);
double rsErfInv(const double x);
double rsOneSampleTTest(const double *data, const unsigned int length, const double mu);
    
void rsRankingResolveTies(double *ranks, const size_t *tiesTrace, const size_t n_ties);
void rsSpearmannRank(double *ranks, const double *data, const size_t n);

double rsComputePValueFromTValue(const double T, const int df);
double rsComputeTValueFromPValue(const double P, const int df);

#ifdef __cplusplus
}
#endif

#endif
