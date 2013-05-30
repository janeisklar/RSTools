/*******************************************************************
 *
 * rsmathutils.cpp
 *
 * 
 * André Hoffmann
 *******************************************************************/

/* C-API that goes into the library and can be accessed without linking anything related to R */
#include "rsmathutils.h"

/* Internal API with all the functions that shouldn't be visible to the outside(anything related to R) */
#include "rsmathutils_internal.h"

void rsLinearRegressionR(int nSamples, double *signal, int nRegressors, double **regressors, double *betas, double *residuals, double *fitted, int verbose)
{
    // create an embedded R instance
    const char * const argv = "";
    static RInside R(0, &argv, true, verbose > 0, false); //argc, argv, loadRcpp, verbose, interactive
    
    // convert inputs to R variables
    R["signal"] = rsMathCreateVector(signal, nSamples);
    
    // build together the regression R command with all regressors
    // > fm <- lm(signal ~ regressor1 + regressor2 + .. + regressorn)
    std::stringstream lmRCommand;
    lmRCommand << "fm <- lm(signal ~ ";
    
    for (int r=0; r<nRegressors; r++) {
        std::stringstream regressorName;
        regressorName << "regressor";
        regressorName << (r+1);
        
        if (regressors != NULL) {
            R[regressorName.str()] = rsMathCreateVector(regressors[r], nSamples);
        }
        
        if (verbose) {
            fprintf(stdout, "Regressor %d: \n", (r+1));
            R.parseEvalQ(regressorName.str());
        }
        
        if ( r>0 ) {
            lmRCommand << " + ";
        }
        
        lmRCommand << regressorName.str();
    }
    lmRCommand << ")";

    if (verbose) {
        fprintf(stdout, "Signal:\n");
        R.parseEvalQ("signal");
    }
    
    // do linear regression
    R.parseEvalQ(lmRCommand.str());
    
    // extract the results and convert them into C variables
    if ( residuals != NULL ) {
        if (verbose) fprintf(stdout, "Residuals:\n");
        Rcpp::NumericVector Rresiduals((SEXP) R.parseEval("residuals(fm)"));
        rsMathCovertVector(Rresiduals, residuals);
    }
    
    if ( fitted != NULL ) {
        if (verbose) fprintf(stdout, "Fitted:\n");
        Rcpp::NumericVector Rfitted((SEXP) R.parseEval("fitted.values(fm)"));
        rsMathCovertVector(Rfitted, fitted);
    }
    if ( betas != NULL ) {
        if (verbose) fprintf(stdout, "Betas:\n");
        Rcpp::NumericVector Rbetas((SEXP) R.parseEval("as.numeric(coefficients(fm))"));
        rsMathCovertVector(Rbetas, betas);
    }
}

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

Rcpp::NumericVector rsMathCreateVector(double *v, int length)
{
    Rcpp::NumericVector V(length);
    for (int i=0; i<length; i++) {
        V(i) = v[i];
    }
    return V;
}

void rsMathCovertVector(Rcpp::NumericVector V, double *target)
{
    for (int i=0; i<V.size(); i++) {
        target[i] = (double)V[i];
    }
}

double* rsMathCovertDoubleVector(Rcpp::NumericVector V)
{
    double *result = (double*)malloc( V.size() * sizeof(double) );
    
    for (int i=0; i<V.size(); i++) {
        result[i] = (double)V[i];
    }
    
    return result;
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