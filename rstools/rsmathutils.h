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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


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

#ifdef USE_OPENBLAS
void openblas_set_num_threads(int num_threads);
#endif
    
#define RSFFTFILTER_CUTOFF 1
#define RSFFTFILTER_SIGMOID 2

#define RSFFTFILTER_ENGINE_GSL 1
#define RSFFTFILTER_ENGINE_FFTW 2
    
#define RSMATRIXCONVERSION_ABSOLUTE 1
#define RSMATRIXCONVERSION_POSITIVE 2
#define RSMATRIXCONVERSION_NEGATIVE 3
#define RSMATRIXCONVERSION_SCALED   4
    
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

struct rsPCAResult {
	gsl_matrix* transformed;
	gsl_vector* eigenvalues;
	gsl_vector* eigenvalues_all;
	gsl_matrix* eigenvectors;
};

struct rsCTPResult {
	gsl_matrix* transformedA;
	gsl_matrix* transformedB;
	gsl_vector* eigenvalues;
	gsl_vector* eigenvalues_all;
	gsl_matrix* eigenvectors;
};

void rsLinearRegression(const int nSamples, const double *signal, const int nRegressors, const double **regressors, double *betas, double *residuals, double *fitted, const int zScoreRegression, const int verbose);
void rsLinearRegressionFilter(const int nSamples, const double *signal, const int nRegressors, const double **regressors, const double sampling_rate, const double f1, const double f2, double *betas, double *residuals, double *fitted, const int verbose);
struct rsFFTFilterParams rsFFTFilterInit(const int T, const long paddedT, const double sampling_rate, const double f1, const double f2, const int rolloff_method, const double rolloff, const BOOL keepMean, const int verbose);
void rsFFTFilter(struct rsFFTFilterParams p, double *data);
void rsFFTFilterFree(struct rsFFTFilterParams p);
void rsFFTSetEngine(int engine);
BOOL rsVoxelInSphere(FloatPoint3D point, FloatPoint3D center, double radius);
BOOL rsVoxelInCube(FloatPoint3D point, FloatPoint3D center, FloatPoint3D dim);
double rsCorrelation(const double* X, const double *Y, const size_t length);
double rsZCorrelation(const double* X, const double* Y, const size_t length);
double rsFastZCorrelation(const double* X, const double* Y, const size_t length);
double rsTCorrelation(const double* X, const double* Y, const size_t length);
double rsMonteCarloZCorrelation(const double* X, const double* Y, const size_t length, const unsigned int repetitions, const unsigned int samplingSize);
double rsDistance(FloatPoint3D A, FloatPoint3D B);
double rsSampleSineWave(const double sampling_rate, const double f, const int t);
double rsSampleCosineWave(const double sampling_rate, const double f, const int t);
double rsSigmoidRolloff(const double nBins, const double rolloff, const double bin);    
double rsSigmoid(const double rolloff, const double x);
long double *rsFirstEigenvector(const double **A, const long n, const long maxIterations, const double precision, const BOOL verbose);
struct rsPCAResult rsGenericPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose);
struct rsPCAResult rsPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose);
struct rsPCAResult rsTPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose);
void rsPCAResultFree(struct rsPCAResult result);
struct rsCTPResult rsCTP(const gsl_matrix* A, const gsl_matrix* B, int nComponents, BOOL verbose);
void rsCTPResultFree(struct rsCTPResult result);
long rsMakePositiveDefiniteSymmetric(gsl_matrix* A);
long double *rsMatrixByVectorProduct(const double **A, const long double *x, const long n, const long m);
long double rsEuclideanNorm(const long double *x, const long n);
void rsScaleVector(long double *x, const long n, const long double factor);
void rsVectorSub(long double *x, const long double *y, const long n);
void rsVectorSwap(long double *x, long double *y, const long n);
long double rsVectorMean(const long double *x, const long n);
BOOL rsVectorContains(const long *x, const long n, const long element);
void rsMatrixConversion(double **A, const long m, const long n, const int mode);

double **d2matrix(int yh, int xh);
int rs_gsl_matrix_fprintf(FILE *stream,gsl_matrix *m,char *fmt);
void rs_gsl_vector_fprintf(FILE *stream, gsl_vector *v, const char* fmt);
void rs_matrix_fprintf(FILE *stream, const double **A, const long m, const long n, const char* fmt);
void rs_vector_fprintfl(FILE *stream, const long double *x, const long n, const char* fmt);
void rs_vector_fprintf(FILE *stream, const double *x, const long n, const char* fmt);
BOOL rsSaveMatrix(const char *filename, const double** A, const long m, const long n);
BOOL rsLoadMatrix(const char *filename, double** A, const long m, const long n);
int rsCountDigits(int n);
    
float  rsFastErfInv(float x);
double rsErfInv(const double x);
double rsOneSampleTTest(const double *data, const unsigned int length, const double mu);
    
void rsRankingResolveTies(double *ranks, const size_t *tiesTrace, const size_t n_ties);
void rsSpearmannRank(double *ranks, const double *data, const size_t n);

gsl_rng *rsGetRandomNumberGenerator();
void rsDestroyRandomNumberGenerator();
    
#ifdef __cplusplus
}
#endif
    
#endif
