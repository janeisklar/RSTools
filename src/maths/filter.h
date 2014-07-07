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

#include "utils.h"
#include "src/nifti/rsniftiutils.h"

#if !defined(__FILTER_MATHUTILS_H)
#define __FILTER_MATHUTILS_H

#if !defined(RS_FFTW_ENABLED)
#define RS_FFTW_ENABLED 0
#endif

#if RS_FFTW_ENABLED == 1
#include <fftw3.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define RSFFTFILTER_CUTOFF 1
#define RSFFTFILTER_SIGMOID 2

#define RSFFTFILTER_ENGINE_GSL 1
#define RSFFTFILTER_ENGINE_FFTW 2

typedef struct {
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
} rsFFTFilterParams;

rsFFTFilterParams* rsFFTFilterInit(const int T, const long paddedT, const double sampling_rate, const double f1, const double f2, const int rolloff_method, const double rolloff, const BOOL keepMean, const int verbose);
void rsFFTFilter(rsFFTFilterParams *p, double *data);
void rsFFTFilterFree(rsFFTFilterParams *p);
void rsFFTSetEngine(int engine);

#ifdef __cplusplus
}
#endif

#endif
