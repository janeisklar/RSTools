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

struct rsPCAResult rsGenericPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose);
struct rsPCAResult rsPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose);
struct rsPCAResult rsTPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose);
void rsPCAResultFree(struct rsPCAResult result);
struct rsCTPResult rsCTP(const gsl_matrix* A, const gsl_matrix* B, int nComponents, BOOL verbose);
void rsCTPResultFree(struct rsCTPResult result);
    
#endif