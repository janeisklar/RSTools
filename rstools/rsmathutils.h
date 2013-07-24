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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include "rsniftiutils.h"

#if !defined(RS_FFTW_ENABLED)
#define RS_FFTW_ENABLED 0
#endif

#if RS_FFTW_ENABLED == 1
#include <fftw3.h>
#endif

#if !defined(__MATHUTILS_H)
#define __MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif
    
#define RSFFTFILTER_CUTOFF 1
#define RSFFTFILTER_SIGMOID 2

#define RSFFTFILTER_ENGINE_GSL 1
#define RSFFTFILTER_ENGINE_FFTW 2
    
struct rsFFTFilterParams {
    int T;
    long paddedT;
    double sampling_rate;
    double f1;
    double f2;
    
    double *frequencyBins;
    double *binAttenuation;
    
    int verbose;
    
    int rolloff_method;
    int rolloff;
    
#if RS_FFTW_ENABLED == 1
    fftw_plan plan_r2hc, plan_hc2r;
#endif
};
    
void rsLinearRegression(const int nSamples, const double *signal, const int nRegressors, const double **regressors, double *betas, double *residuals, double *fitted, const int verbose);
void rsLinearRegressionFilter(const int nSamples, const double *signal, const int nRegressors, const double **regressors, const double sampling_rate, const double f1, const double f2, double *betas, double *residuals, double *fitted, const int verbose);
struct rsFFTFilterParams rsFFTFilterInit(const int T, const long paddedT, const double sampling_rate, const double f1, const double f2, const int rolloff_method, const double rolloff, const int verbose);
void rsFFTFilter(struct rsFFTFilterParams p, double *data);
void rsFFTFilterFree(struct rsFFTFilterParams p);
void rsFFTSetEngine(int engine);
BOOL rsVoxelInSphere(FloatPoint3D point, FloatPoint3D center, double radius);
BOOL rsVoxelInCube(FloatPoint3D point, FloatPoint3D center, FloatPoint3D dim);
double rsCorrelation(const double* X, const double *Y, const size_t length);
double rsZCorrelation(const double* X, const double* Y, const size_t length);
double rsDistance(FloatPoint3D A, FloatPoint3D B);
double rsSampleSineWave(const double sampling_rate, const double f, const int t);
double rsSampleCosineWave(const double sampling_rate, const double f, const int t);
double rsSigmoidRolloff(const double nBins, const double rolloff, const double bin);    
double rsSigmoid(const double rolloff, const double x);
gsl_vector *rsFirstEigenvector(gsl_matrix* A, long maxIterations, double precision);
double **d2matrix(int yh, int xh);
    
#ifdef __cplusplus
}
#endif
    
#endif