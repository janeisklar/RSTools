#include "regression.h"

rsMultifitLinearWorkspace* rsPrepareLinearRegression(const int nSamples, const int nRegressors, const double **regressors, const int zScoreRegression)
{
    // create workspace
    rsMultifitLinearWorkspace *work = rsMultifitLinearWorkspaceAlloc(nSamples, nRegressors);
    
    // create regressor matrix / design matrix
    gsl_matrix *X   = work->X; 
    
    for (int r=0; r<nRegressors; r++) {
        for (int t=0; t<nSamples; t++) {
            gsl_matrix_set(X,t,r,regressors[r][t]);
        }
    }
    
    // remove the regressors' mean
    for (int i=0; i<nRegressors; i=i+1) {
            
        gsl_vector_view regressor = gsl_matrix_column(X, i);
        
        const double meanX = gsl_stats_mean(regressor.vector.data, 1, nSamples);
        
        gsl_vector_add_constant(&regressor.vector, -1.0 * meanX); // remove mean

        if ( zScoreRegression ) {
            const double sX = gsl_stats_sd_m(regressor.vector.data, 1, nSamples, meanX);
            gsl_vector_scale(&regressor.vector, 1 / sX); // divide by st.dev
        }
    }
    
    // precompute the pseudo-inverse of the regressor matrix
    rsPrepareMultifitLinearSVD(X, 1, work);
    
    return work;
}

void rsLinearRegression(rsMultifitLinearWorkspace *work, const double *signal, double *betas, double *residuals, double *fitted, const int zScoreRegression, const int verbose)
{
    const int nSamples    = work->X->size1;
    const int nRegressors = work->X->size2;
    
    // convert inputs to gsl variables
    gsl_vector *y   = gsl_vector_alloc(nSamples);                 // input/signal
    gsl_vector *b   = gsl_vector_alloc(nRegressors);              // betas
    gsl_matrix *cov = gsl_matrix_alloc(nRegressors, nRegressors); // covariance
    gsl_vector *res = gsl_vector_alloc(nSamples);                 // residuals
    
    // fill them
    for (int i=0; i < nSamples; i++){
        gsl_vector_set(y,i,signal[i]);
    }
    
    // transform to z-scores if requested
    // according to http://www.personal.kent.edu/~jortiz/earthstats/topic06notes.html#defining
    if ( zScoreRegression ) {
        
        // transform the input signal
        const double meanY = gsl_stats_mean(signal, 1, nSamples);
        const double sY = gsl_stats_sd_m(signal, 1, nSamples, meanY);
        
        gsl_vector_add_constant(y, -1.0 * meanY); // remove mean
        gsl_vector_scale(y, 1 / sY); // divide by st.dev
    }
    
    // execute linear regression
    double chisq;
    size_t rank;
    int success = rsMultifitLinearSVD (y, &rank, b, cov, &chisq, work);
    
    // compute residuals
    gsl_multifit_linear_residuals(work->X, y, b, res);
    
    // compute fits
    if ( fitted != NULL ) {
        
        double err;
        
        for (int t=0; t<nSamples; t=t+1) {
            gsl_vector_const_view y = gsl_matrix_const_row(work->X, t);
            gsl_multifit_linear_est(&y.vector, b, cov, &fitted[t], &err);
        }
    }
    
    // convert back to basic C variables
    for (int t=0; t < nSamples; t++) {
        residuals[t] = gsl_vector_get(res, t);
    }
    
    for (int r=0; r < nRegressors; r++) {
        betas[r] = gsl_vector_get(b, r);
    }
    
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(b);
    gsl_vector_free(res);
}

/* Fit
 *
 *  y = X c
 *
 *  where X is an M x N matrix of M observations for N variables.
 * 
 * This is a modified copy of gsl's multifit_linear_svd
 */
int rsPrepareMultifitLinearSVD (const gsl_matrix * X, int balance, rsMultifitLinearWorkspace * work) {
    
    gsl_matrix *A   = work->w->A;
    gsl_matrix *Q   = work->w->Q;
    gsl_matrix *QSI = work->w->QSI;
    gsl_vector *S   = work->w->S;
    gsl_vector *xt  = work->w->xt;
    gsl_vector *D   = work->w->D;

    /* Copy X to workspace,  A <= X */
    gsl_matrix_memcpy(A, X);
    
    /* Copy X to workspace,  A <= X */
    gsl_matrix_memcpy(work->X, X);

    /* Balance the columns of the matrix A if requested */
    if (balance) {
      gsl_linalg_balance_columns(A, D);
    } else {
      gsl_vector_set_all(D, 1.0);
    }

    /* Decompose A into U S Q^T */
    gsl_linalg_SV_decomp_mod(A, QSI, Q, S, xt);

    return GSL_SUCCESS;
}

/* Fit
 *
 *  y = X c
 *
 *  where X is an M x N matrix of M observations for N variables.
 * 
 * This is a modified copy of gsl's multifit_linear_svd
 */
int rsMultifitLinearSVD (const gsl_vector * y,
                         size_t * rank,
                         gsl_vector * c,
                         gsl_matrix * cov,
                         double *chisq, 
                         rsMultifitLinearWorkspace * work)
{

    const double tol = GSL_DBL_EPSILON;
    
    if (work->w->A->size1 != y->size) {
        GSL_ERROR("number of observations in y does not match rows of matrix X", GSL_EBADLEN);
    } else if (work->w->A->size2 != c->size) {
        GSL_ERROR ("number of parameters c does not match columns of matrix X", GSL_EBADLEN);
    } else if (cov->size1 != cov->size2) {
        GSL_ERROR ("covariance matrix is not square", GSL_ENOTSQR);
    } else if (c->size != cov->size1) {
        GSL_ERROR("number of parameters does not match size of covariance matrix", GSL_EBADLEN);
    }
    #if GSL_MAJOR_VERSION < 2
    else if (work->w->A->size1 != work->w->n || work->w->A->size2 != work->w->p) {
    #else
    else if (work->w->A->size1 != work->w->nmax || work->w->A->size2 != work->w->pax) {
    #endif
        GSL_ERROR("size of workspace does not match size of observation matrix", GSL_EBADLEN);
    } else {
        size_t i, j, p_eff;

        gsl_matrix *A   = work->w->A;
        gsl_matrix *Q   = work->w->Q;
        gsl_matrix *QSI = work->w->QSI;
        gsl_vector *S   = work->w->S;
        gsl_vector *xt  = work->w->xt;
        gsl_vector *D   = work->w->D;
        gsl_matrix *X   = work->X;
        
        const size_t n = A->size1;
        const size_t p = A->size2;
        
        /* Solve y = A c for c */
        
        gsl_blas_dgemv(CblasTrans, 1.0, A, y, 0.0, xt);
        
        /* Scale the matrix Q,  Q' = Q S^-1 */
        
        gsl_matrix_memcpy(QSI, Q);
        
        {
            double alpha0 = gsl_vector_get(S, 0);
            p_eff = 0;
            
            for (j = 0; j < p; j++) {
                gsl_vector_view column = gsl_matrix_column(QSI, j);
                double alpha = gsl_vector_get(S, j);
            
                if (alpha <= tol * alpha0) {
                    alpha = 0.0;
                } else {
                    alpha = 1.0 / alpha;
                    p_eff++;
                }
            
                gsl_vector_scale(&column.vector, alpha);
              }
            
            *rank = p_eff;
        }
        
        gsl_vector_set_zero(c);
        
        gsl_blas_dgemv(CblasNoTrans, 1.0, QSI, xt, 0.0, c);
        
        /* Unscale the balancing factors */
        
        gsl_vector_div(c, D);
        
        /* Compute chisq, from residual r = y - X c */
        
        {
          double s2 = 0, r2 = 0;
        
          for (i = 0; i < n; i++) {
              double yi = gsl_vector_get(y, i);
              gsl_vector_const_view row = gsl_matrix_const_row(X, i);
              double y_est, ri;
              gsl_blas_ddot(&row.vector, c, &y_est);
              ri = yi - y_est;
              r2 += ri * ri;
          }
        
          s2 = r2 / (n - p_eff);   /* p_eff == rank */
        
          *chisq = r2;
        
          /* Form variance-covariance matrix cov = s2 * (Q S^-1) (Q S^-1)^T */
        
          for (i = 0; i < p; i++) {
              gsl_vector_view row_i = gsl_matrix_row(QSI, i);
              double d_i = gsl_vector_get(D, i);
        
              for (j = i; j < p; j++) {
                  gsl_vector_view row_j = gsl_matrix_row(QSI, j);
                  double d_j = gsl_vector_get(D, j);
                  double s;
        
                  gsl_blas_ddot(&row_i.vector, &row_j.vector, &s);
        
                  gsl_matrix_set(cov, i, j, s * s2 /(d_i * d_j));
                  gsl_matrix_set (cov, j, i, s * s2 / (d_i * d_j));
                }
            }
        }
        
        return GSL_SUCCESS;
    }
}


rsMultifitLinearWorkspace* rsMultifitLinearWorkspaceAlloc(size_t n, size_t p) {
    
    rsMultifitLinearWorkspace *w = (rsMultifitLinearWorkspace *) rsMalloc(sizeof(rsMultifitLinearWorkspace));
    gsl_multifit_linear_workspace *gw = gsl_multifit_linear_alloc(n,p);
    
    w->w = gw;
    w->X = gsl_matrix_alloc(n, p);
    
    return w;
}

void rsMultifitLinearWorkspaceFree(rsMultifitLinearWorkspace* work) {
    gsl_matrix_free(work->X);
    gsl_multifit_linear_free(work->w);
    rsFree(work);
}

rsMultifitLinearWorkspace* rsMultifitLinearWorkspaceClone(rsMultifitLinearWorkspace* work) {
    rsMultifitLinearWorkspace *w = rsMultifitLinearWorkspaceAlloc(work->X->size1, work->X->size2);
    
    gsl_matrix * A;
    gsl_matrix * Q;
    gsl_matrix * QSI;
    gsl_vector * S;
    gsl_vector * t;
    gsl_vector * xt;
    gsl_vector * D;
    
    gsl_matrix_memcpy(w->X,      work->X);
    gsl_matrix_memcpy(w->w->A,   work->w->A);
    gsl_matrix_memcpy(w->w->Q,   work->w->Q);
    gsl_matrix_memcpy(w->w->QSI, work->w->QSI);
    gsl_vector_memcpy(w->w->S,   work->w->S);
    gsl_vector_memcpy(w->w->t,   work->w->t);
    gsl_vector_memcpy(w->w->xt,  work->w->xt);
    gsl_vector_memcpy(w->w->D,   work->w->D);
    
    return w;
}
