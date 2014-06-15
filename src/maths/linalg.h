#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "utils.h"

#if !defined(__LINALG_MATHUTILS_H)
#define __LINALG_MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif 
    
#define RSMATRIXCONVERSION_ABSOLUTE 1
#define RSMATRIXCONVERSION_POSITIVE 2
#define RSMATRIXCONVERSION_NEGATIVE 3
#define RSMATRIXCONVERSION_SCALED   4

long double *rsFirstEigenvector(const double **A, const long n, const long maxIterations, const double precision, const BOOL verbose);
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

#ifdef __cplusplus
}
#endif

#endif