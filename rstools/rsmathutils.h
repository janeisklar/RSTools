#include <stdio.h>
#include <strings.h>

#if !defined(__MATHUTILS_H)
#define __MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

void rsLinearRegression(int nSamples, double *signal, int nRegressors, double **regressors, double *betas, double *residuals, double *fitted, int verbose);
double **d2matrix(int yh, int xh);

#ifdef __cplusplus
}
#endif
    
#endif