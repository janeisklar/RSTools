#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "nifti/rsniftiutils.h"

#if !defined(__UTILS_MATHUTILS_H)
#define __UTILS_MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif    

double rsSampleSineWave(const double sampling_rate, const double f, const int t);
double rsSampleCosineWave(const double sampling_rate, const double f, const int t);
double rsSigmoidRolloff(const double nBins, const double rolloff, const double bin);    
double rsSigmoid(const double rolloff, const double x);
int rsCountDigits(int n);

gsl_rng *rsGetRandomNumberGenerator();
void rsDestroyRandomNumberGenerator();

double **rsLoadRegressors(char *path, long *nRegressors, long *nValues, double constantFactor);
BOOL rsReadline(FILE *f, char *line, int *length);
double *rsParseRegressorLine(char *line, long *nRegressors);

#ifdef __cplusplus
}
#endif

#endif
