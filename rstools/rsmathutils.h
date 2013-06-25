#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_statistics_double.h>
#include "rsniftiutils.h"

#if !defined(__MATHUTILS_H)
#define __MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

void rsLinearRegression(const int nSamples, const double *signal, const int nRegressors, const double **regressors, double *betas, double *residuals, double *fitted, const int verbose);
void rsLinearRegressionFilter(const int nSamples, const double *signal, const int nRegressors, const double **regressors, const double sampling_rate, const double f1, const double f2, double *betas, double *residuals, double *fitted, const int verbose);
void rsFFTFilter(double *data, const int T, const double sampling_rate, const double f1, const double f2, const int verbose);
BOOL rsVoxelInSphere(FloatPoint3D point, FloatPoint3D center, double radius);
BOOL rsVoxelInCube(FloatPoint3D point, FloatPoint3D center, FloatPoint3D dim);
double rsCorrelation(const double* X, const double *Y, const size_t length);
double rsZCorrelation(const double* X, const double* Y, const size_t length);
double rsDistance(FloatPoint3D A, FloatPoint3D B);
double rsSampleSineWave(const double sampling_rate, const double f, const int t);
double rsSampleCosineWave(const double sampling_rate, const double f, const int t);
double **d2matrix(int yh, int xh);

#ifdef __cplusplus
}
#endif
    
#endif