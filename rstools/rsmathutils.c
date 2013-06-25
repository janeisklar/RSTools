/*******************************************************************
 *
 * rsmathutils.c
 *
 * 
 * André Hoffmann
 *******************************************************************/

#include "rsmathutils.h"

void rsLinearRegression(const int nSamples, const double *signal, const int nRegressors, const double **regressors, double *betas, double *residuals, double *fitted, const int verbose)
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
    regressors2 = d2matrix(nRegressors2, nSamples);
    
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
        verbose
    );
    
    // Copy only the beta values of the original regressors
    for (int i=0; i<nRegressors; i=i+1) {
        betas[i] = betas2[i];
    }
    
    free(regressors2[0]);
    free(regressors2);
}

void rsFFTFilter(double *data, const int T, const double sampling_rate, const double f1, const double f2, const int verbose) {
    
    /* Compute the frequency of the spectral bins */
    double F[T];
    for (int i=0; i<T; i=i+1) {
        F[i] = 0.0;
        
        if ( i > 0 ) {
            i      = i+1; // skip the bin holding the real part
            F[i]   = i/(2 * T * sampling_rate); // set the frequency of the complex part's bin
            F[i-1] = F[i]; // copy it to the real part's bin
        }        
    }
    
    /* Compute which bins to keep */
    int i1=-1, i2=-1;
    for (int i=1; i<=T; i=i+2) {
        if (F[i] > f1) {
            i1 = (i-2 < 0) ? 0 : i-2;
            break;
        }
    }
    for (int i=1; i<=T; i=i+2) {
        if (F[i] > f2) {
            i2 = (i+1 >= T) ? T-1 : i+1;
            break;
        }
    }
    if ( i1 < 0 ) i1 = (T % 2==0) ? T-1 : T-2; /* in the even case the complex part */
    if ( i2 < 0 ) i2 = (T % 2==0) ? T-1 : T-2; /* of the last bin is not stored     */
    
    if (verbose) printf("Bandpass range: %.4fHz(%d)..%.4fHz(%d)\n", F[i1], i1, F[i2], i2);
    
    /* Prepare FFT Filtering */
    gsl_fft_real_wavetable        *real;
    gsl_fft_halfcomplex_wavetable *hc;
    gsl_fft_real_workspace        *work;
    
    work = gsl_fft_real_workspace_alloc(T);
    real = gsl_fft_real_wavetable_alloc(T);
    hc   = gsl_fft_halfcomplex_wavetable_alloc(T);
    
    /* FFT */
    gsl_fft_real_transform(data, 1, T, real, work);
    
    /* Remove undesired frequency bins */
    for (int i = 0; i<T; i=i+1) {
        if ( i < i1 || i > i2 ) {
            data[i] = 0;
        }
    }
    
    /* Inverse FFT */
    gsl_fft_halfcomplex_inverse(data, 1, T, hc, work);
    
    /* Free memory */
    gsl_fft_real_wavetable_free(real);
    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);
}

double rsZCorrelation(const double* X, const double* Y, const size_t length)
{
    /* 
     * This implementation is based on:
     * Weissenbacher, Andreas, et al. 
     * "Correlations and anticorrelations in resting-state functional connectivity MRI: a quantitative comparison of preprocessing strategies."
     * Neuroimage 47.4 (2009): 1408-1416.
     */
    
    const long double N = length;
    const long double Nm1 = length - 1;
    
    // compute means
    long double meanX = 0.0;
    long double meanY = 0.0;
    
    for (unsigned int i=0; i<length; i=i+1) {
        meanX = meanX + (long double)X[i] / N;
        meanY = meanY + (long double)Y[i] / N;
    }
    
    // compute standard scores(std. dev)
    long double sX = 0.0;
    long double sY = 0.0;

    for (unsigned int i=0; i<length; i=i+1) {
        sX = sX + powl((long double)X[i]-meanX, 2.0) / Nm1;
        sY = sY + powl((long double)Y[i]-meanY, 2.0) / Nm1;
    }
    
    sX = sqrtl(sX);
    sY = sqrtl(sY);
    
    // compute correlation coeeficient
    long double r = 0.0;
    const long double norm = sX * sY * Nm1;
    
    for (unsigned int i=0; i<length; i=i+1) {
        r = r + (long double)((long double)X[i]-meanX) * ((long double)Y[i]-meanY) / norm;
    }
    
    // Fisher's r-to-z transformation
    const long double one = 1.0;
    const long double half = 0.5;
    return half * logl( (one + r) / (one - r) );
}

double rsCorrelation(const double* X, const double* Y, const size_t length)
{
    return gsl_stats_correlation(X, 1,
                                 Y, 1,
                                 length);
}

double rsDistance(FloatPoint3D A, FloatPoint3D B)
{
    return pow(
        (
            pow(A.x-B.x,2.0) +
            pow(A.y-B.y,2.0) +
            pow(A.z-B.z,2.0)
        ),
        1.0/2.0
    );
}

BOOL rsVoxelInSphere(FloatPoint3D point, FloatPoint3D center, double radius)
{
    return rsDistance(point, center) <= radius;
}

BOOL rsVoxelInCube(FloatPoint3D point, FloatPoint3D center, FloatPoint3D dim)
{
    return fabs(point.x-center.x) <= (dim.x / 2.0) &&
           fabs(point.y-center.y) <= (dim.y / 2.0) &&
           fabs(point.z-center.z) <= (dim.z / 2.0);
}

double rsSampleSineWave(const double sampling_rate, const double f, const int t)
{
    return sin((double)t * 2.0 * M_PI * f * sampling_rate);
}

double rsSampleCosineWave(const double sampling_rate, const double f, const int t)
{
    return cos((double)t * 2.0 * M_PI * f * sampling_rate);
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