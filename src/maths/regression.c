#include "regression.h"

void rsLinearRegression(const int nSamples, const double *signal, const int nRegressors, const double **regressors, double *betas, double *residuals, double *fitted, const int zScoreRegression, const int verbose)
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
    
    // transform to z-scores if requested
    // according to http://www.personal.kent.edu/~jortiz/earthstats/topic06notes.html#defining
    if ( zScoreRegression ) {
        
        // transform the input signal
        const double meanY = gsl_stats_mean(signal, 1, nSamples);
        const double sY = gsl_stats_sd_m(signal, 1, nSamples, meanY);
        
        gsl_vector_add_constant(y, -1.0 * meanY); // remove mean
        gsl_vector_scale(y, 1 / sY); // divide by st.dev
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
    
    // execute linear regression
    double chisq;
    int success = gsl_multifit_linear(X, y, b, cov, &chisq, work);
    
    // compute residuals
    gsl_multifit_linear_residuals(X, y, b, res);
    
    // compute fits
    if ( fitted != NULL ) {
        
        double err;
        
        for (int t=0; t<nSamples; t=t+1) {
            gsl_vector_const_view y = gsl_matrix_const_row(X, t);
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
    
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(b);
    gsl_vector_free(res);
    gsl_multifit_linear_free(work);
}

void rsLinearRegressionFilter(
    const int nSamples,
    const double *signal,
    const int nRegressors,
    const double **regressors,
    const double sampling_rate,
    const double f1,
    const double f2,
    double *betas,
    double *residuals,
    double *fitted,
    const int verbose) {

    const double nyquist_frequency = 1.0 / (2.0 * sampling_rate);
    const double bin_width = nyquist_frequency / (nSamples/3);
    
    // Compute the number of frequency regressors that will be added to the existing regressors
    const int nFrequencyBinsLow  = (int)floor(f1 / bin_width - 1.0);                       // number of frequency bins before the lower cutoff frequency
    const int nFrequencyBinsHigh = (int)floor((nyquist_frequency - f2) / bin_width - 1.0); // number of frequency bins after the highter cutoff frequency
    const int nFrequencyBins = nFrequencyBinsLow + nFrequencyBinsHigh;                     // number of frequency bins in total
    const int nFrequencyRegressors = nFrequencyBins * 2;                                   // number of frequency regressors (both sine and cosine)
    
    // Compute the frequencies for the bins
    double frequencyBins[nFrequencyBins];
    
    for ( int i=0; i<nFrequencyBins; i=i+1 ) {
        if ( (i < nFrequencyBinsLow) ) {
            frequencyBins[i] = (i + 1) * bin_width;
        } else {
            frequencyBins[i] = (i - nFrequencyBinsLow + 2) * bin_width;
        }
    }
    
    // Create new regressor matrix
    const int nRegressors2 = nRegressors + nFrequencyRegressors;
    double **regressors2;
    regressors2 = d2matrix(nRegressors2-1, nSamples-1);
    
    for (int i=0; i<nRegressors2; i=i+1) {
        
        if ( i < nRegressors ) {
            for (int t=0; t<nSamples; t=t+1) {
                regressors2[i][t] = regressors[i][t];
            }
        } else if ( i < (nRegressors + nFrequencyBins) ) {
            const int j = i - nRegressors;
            for (int t=0; t<nSamples; t=t+1) {
                regressors2[i][t] = rsSampleSineWave(sampling_rate, frequencyBins[j], t);
            }
        } else {
            const int j = i - nRegressors - nFrequencyBins;
            for (int t=0; t<nSamples; t=t+1) {
                regressors2[i][t] = rsSampleCosineWave(sampling_rate, frequencyBins[j], t);
            }
        }
    }
    
    double betas2[nRegressors2];
    
    // Regress
    rsLinearRegression(
        nSamples,
        signal,
        nRegressors2,
        (const double**)regressors2,
        betas2,
        residuals,
        fitted,
        0,
        verbose
    );
    
    // Copy only the beta values of the original regressors
    for (int i=0; i<nRegressors; i=i+1) {
        betas[i] = betas2[i];
    }
    
    free(regressors2[0]);
    free(regressors2);
}
