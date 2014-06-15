#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "linalg.h"
#include "utils.h"

#if !defined(__MULTIVAR_ANALYSIS_MATHUTILS_H)
#define __MULTIVAR_ANALYSIS_MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	gsl_matrix* transformed;
	gsl_vector* eigenvalues;
	gsl_vector* eigenvalues_all;
	gsl_matrix* eigenvectors;
} rsPCAResult;

typedef struct {
	gsl_matrix* transformedA;
	gsl_matrix* transformedB;
	gsl_vector* eigenvalues;
	gsl_vector* eigenvalues_all;
	gsl_matrix* eigenvectors;
} rsCSPResult;

rsPCAResult *rsGenericPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose);
rsPCAResult *rsPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose);
rsPCAResult *rsTPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose);
void rsPCAResultFree(rsPCAResult *result);
rsCSPResult *rsCSP(const gsl_matrix* A, const gsl_matrix* B, int nComponents, BOOL verbose);
void rsCSPResultFree(rsCSPResult *result);

#ifdef __cplusplus
}
#endif

#endif
