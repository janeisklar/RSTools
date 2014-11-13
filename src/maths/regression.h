#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "linalg.h"

#if !defined(__REGRESSION_MATHUTILS_H)
#define __REGRESSION_MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif 
    

typedef struct  {
    gsl_multifit_linear_workspace* w;
    gsl_matrix* X;
} rsMultifitLinearWorkspace;

rsMultifitLinearWorkspace* rsPrepareLinearRegression(const int nSamples, const int nRegressors, const double **regressors, const int zScoreRegression);
void rsLinearRegression(rsMultifitLinearWorkspace *work, const double *signal, double *betas, double *residuals, double *fitted, const int zScoreRegression, const int verbose);

int rsPrepareMultifitLinearSVD (const gsl_matrix * X, int balance, rsMultifitLinearWorkspace * work);
int rsMultifitLinearSVD (const gsl_vector * y,
                         size_t * rank,
                         gsl_vector * c,
                         gsl_matrix * cov,
                         double *chisq, 
                         rsMultifitLinearWorkspace * work);

rsMultifitLinearWorkspace* rsMultifitLinearWorkspaceAlloc(size_t n, size_t p);
void rsMultifitLinearWorkspaceFree(rsMultifitLinearWorkspace* work);
rsMultifitLinearWorkspace* rsMultifitLinearWorkspaceClone(rsMultifitLinearWorkspace* work);

#ifdef __cplusplus
}
#endif

#endif
