/*******************************************************************
 *
 * rsmathutils.cpp
 *
 * 
 * André Hoffmann
 *******************************************************************/

#include "rsmathutils.h"

void rsLinearRegression(int nSamples, double *signal, int nRegressors, double **regressors, double *betas, double *residuals, double *fitted, int verbose)
{
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(nSamples, nRegressors);
    
    // convert inputs to gsl variables
    gsl_vector *y   = gsl_vector_alloc(nSamples);                 // input/signal
    gsl_matrix *X   = gsl_matrix_alloc(nSamples, nRegressors);    // regressors
    gsl_vector *b   = gsl_vector_alloc(nRegressors);              // betas
    gsl_matrix *cov = gsl_matrix_alloc(nRegressors, nRegressors); // covariance
    gsl_vector *res = gsl_vector_alloc(nSamples);                 // residuals
    
    // fill them
    for (int i=0; i < nSamples; i++){
        gsl_vector_set(y,i,signal[i]);
    }
    
    for (int r=0; r<nRegressors; r++) {
        for (int t=0; t<nSamples; t++) {
            gsl_matrix_set(X,t,r,regressors[r][t]);
        }
    }
    
    // execute linear regression
    double chisq;
    int success = gsl_multifit_linear(X, y, b, cov, &chisq, work);
    
    // compute residuals
    gsl_multifit_linear_residuals(X, y, b, res);
    
    // convert back to basic C variables
    for (int t=0; t < nSamples; t++) {
        residuals[t] = gsl_vector_get(res, t);
    }
    
    for (int r=0; r < nRegressors; r++) {
        betas[r] = gsl_vector_get(b, r);
    }
    
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(b);
    gsl_vector_free(res);
    gsl_multifit_linear_free(work);
}

double **d2matrix(int yh, int xh)
{
    int j;
    int nrow = yh+1;
    int ncol = xh+1;
    double **t;
    
    
    /** allocate pointers to ydim */
    t=(double **) malloc((size_t)((nrow)*sizeof(double*)));
    
    /** allocate pointers for zdim */
    t[0]=(double *) malloc((size_t)((ncol*nrow)*sizeof(double)));
    
    /** point everything to the data blob */
    for(j=1;j<nrow;j++) t[j]=t[j-1]+ncol;
    
    return t;
}