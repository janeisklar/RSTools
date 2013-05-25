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

void rsLinearRegression(int nSamples, double *signal, int nRegressors, double **regressors, double *betas, double *residuals, double *fitted, int verbose)
{
    // create an embedded R instance
    const char * const argv = "";
    RInside R(0, &argv, true, verbose > 0, false); //argc, argv, loadRcpp, verbose, interactive
    
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
        
        R[regressorName.str()] = rsMathCreateVector(regressors[r], nSamples);
        
        if ( r>0 ) {
            lmRCommand << " + ";
        }
        
        lmRCommand << regressorName.str();
    }
    lmRCommand << ")";
           
    // do linear regression
    R.parseEvalQ(lmRCommand.str());
    
    // extract the results
    Rcpp::NumericVector Rresiduals( (SEXP) R.parseEval("residuals(fm)"));
    Rcpp::NumericVector Rfitted(    (SEXP) R.parseEval("fitted.values(fm)"));
    Rcpp::NumericVector Rbetas(     (SEXP) R.parseEval("as.numeric(coefficients(fm))"));
    
    // convert the results into C variables
    rsMathCovertVector(Rresiduals, residuals);
    rsMathCovertVector(Rfitted,    fitted);
    rsMathCovertVector(Rbetas,     betas);
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