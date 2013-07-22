/*******************************************************************
 *
 * rsmathutils.c
 *
 * 
 * André Hoffmann
 *******************************************************************/

#include "rsmathutils.h"

static rsFFTFilterEngine = RSFFTFILTER_ENGINE_GSL;

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

struct rsFFTFilterParams rsFFTFilterInit(const int T, const long paddedT, const double sampling_rate, const double f1, const double f2, const int rolloff_method, const double rolloff, const int verbose) {
    
    struct rsFFTFilterParams p;
    double *F = malloc(paddedT * sizeof(double));
    double *attenuation = malloc(paddedT*sizeof(double));

#if RS_FFTW_ENABLED == 1
    if (rsFFTFilterEngine == RSFFTFILTER_ENGINE_FFTW) {
    
        /* Compute the frequency of the spectral bins */
        F[0] = 0.0;
        
        // set the frequency of the real part's bin
        for (int i=1; i<=(int)(paddedT/2); i=i+1) {
            F[i]   = i/(2 * paddedT * sampling_rate);
        }
        
        // set the frequency of the complex part's bin
        for (int i=1; i<=(int)((paddedT+1)/2-1); i=i+1) {
            F[paddedT-i] = i/(2 * paddedT * sampling_rate);
        }
        
        /* Compute which bins to keep */
        
        // real part left
        int ri1=-1;
        for (int i=1; i<=(int)(paddedT/2); i=i+2) {
            if (F[i] > f1) {
                ri1 = (i-1 < 0) ? 0 : i-1;
                break;
            }
        }
        
        // real part right
        int ri2=-1;
        for (int i=1; i<=(int)(paddedT/2); i=i+2) {
            if (F[i] > f2) {
                ri2 = (i-1 < 0) ? 0 : i-1;
                break;
            }
        }
        
        // complex part left
        int ci1=-1;
        for (int i=1; i<=(int)((paddedT+1)/2-1); i=i+2) {
            int index = paddedT-i;
            if (F[index] > f1) {
                ci1 = (index+1 > paddedT-1) ? paddedT-1 : index+1;
                break;
            }
        }
        
        // complex part right
        int ci2=-1;
        for (int i=1; i<=(int)((paddedT+1)/2-1); i=i+2) {
            int index = paddedT-i;
            if (F[index] > f2) {
                ci2 = (index+1 > paddedT-1) ? paddedT-1 : index+1;
                break;
            }
        }
        
        if (verbose) printf("Bandpass range: real %.4fHz(%d)..%.4fHz(%d) / complex %.4fHz(%d)..%.4fHz(%d) \n", F[ri1], ri1, F[ri2], ri2, F[ci1], ci1, F[ci2], ci2);
        
        /* Init attenuation of the bins */
        for (int i = 0; i<=(int)(paddedT/2); i=i+1) {
            attenuation[i] = (i < ri1 || i > ri2) ? 0.0 : 1.0;
        }
        for (int i = (int)(paddedT/2)+1; i<paddedT; i=i+1) {
            attenuation[i] = (i > ci1 || i < ci2) ? 0.0 : 1.0;
            //attenuation[i] = 1.0;
        }

        /*
        if ( rolloff_method == RSFFTFILTER_SIGMOID ) {
            // Add sigmoid attenuation to the lower range
            attenuation[0] = rsSigmoidRolloff(paddedT, rolloff, -i1);
            for (int i = 2; i<i1; i=i+2) {
                attenuation[i] = rsSigmoidRolloff(paddedT, rolloff, i-i1);
                attenuation[i-1] = attenuation[i];
            }
            
            // Add sigmoid attenuation to the upper range
            for (int i = i2+2; i<paddedT; i=i+2) {
                attenuation[i]   = rsSigmoidRolloff(paddedT, rolloff, i-i2);
                attenuation[i-1] = attenuation[i];
            }
            
            if ( paddedT % 2==0 ) {
                attenuation[paddedT-1] = rsSigmoidRolloff(paddedT, rolloff, paddedT-i2-1);
            }
        }
        */
        
        /* Create FFTW transformation plan */
        double *in  = (double*) fftw_malloc(sizeof(double) * paddedT);
        double *out = (double*) fftw_malloc(sizeof(double) * paddedT);
        
        p.plan_r2hc = fftw_plan_r2r_1d(paddedT, in, out, FFTW_R2HC, FFTW_PATIENT);
        p.plan_hc2r = fftw_plan_r2r_1d(paddedT, out, in, FFTW_HC2R, FFTW_PATIENT);

        fftw_free(in);
        fftw_free(out);
    } else if (rsFFTFilterEngine == RSFFTFILTER_ENGINE_GSL) {
#endif
        /* Compute the frequency of the spectral bins */
        for (int i=0; i<paddedT; i=i+1) {
            F[i] = 0.0;
            
            if ( i > 0 ) {
                i      = i+1; // skip the bin holding the real part
                F[i]   = i/(2 * paddedT * sampling_rate); // set the frequency of the complex part's bin
                F[i-1] = F[i]; // copy it to the real part's bin
            }
        }
        
        /* Compute which bins to keep */
        int i1=-1, i2=-1;
        for (int i=1; i<=paddedT; i=i+2) {
            if (F[i] > f1) {
                i1 = (i-2 < 0) ? 0 : i-2;
                break;
            }
        }
        for (int i=1; i<=paddedT; i=i+2) {
            if (F[i] > f2) {
                i2 = (i+1 >= paddedT) ? paddedT-1 : i+1;
                break;
            }
        }
        if ( i1 < 0 ) i1 = (paddedT % 2==0) ? paddedT-1 : paddedT-2; /* in the even case the complex part.. */
        if ( i2 < 0 ) i2 = (paddedT % 2==0) ? paddedT-1 : paddedT-2; /* ..of the last bin is not stored     */
        
        if (verbose) printf("Bandpass range: %.4fHz(%d)..%.4fHz(%d)\n", F[i1], i1, F[i2], i2);
        
        /* Init attenuation of the bins */
        for (int i = 0; i<paddedT; i=i+1) {
            attenuation[i] = (i < i1 || i > i2) ? 0.0 : 1.0;
        }
        
        if ( rolloff_method == RSFFTFILTER_SIGMOID ) {
            // Add sigmoid attenuation to the lower range
            attenuation[0] = rsSigmoidRolloff(paddedT, rolloff, -i1);
            for (int i = 2; i<i1; i=i+2) {
                attenuation[i] = rsSigmoidRolloff(paddedT, rolloff, i-i1);
                attenuation[i-1] = attenuation[i];
            }
            
            // Add sigmoid attenuation to the upper range
            for (int i = i2+2; i<paddedT; i=i+2) {
                attenuation[i]   = rsSigmoidRolloff(paddedT, rolloff, i-i2);
                attenuation[i-1] = attenuation[i];
            }
            
            if ( paddedT % 2==0 ) {
                attenuation[paddedT-1] = rsSigmoidRolloff(paddedT, rolloff, paddedT-i2-1);
            }
        }
#if RS_FFTW_ENABLED == 1
    }
#endif
    
    p.frequencyBins  = F;
    p.binAttenuation = attenuation;
    p.f1             = f1;
    p.f2             = f2;
    p.verbose        = verbose;
    p.sampling_rate  = sampling_rate;
    p.T              = T;
    p.paddedT        = paddedT;
    p.rolloff_method = rolloff_method;
    p.rolloff        = rolloff;
    
    return p;
}

#if RS_FFTW_ENABLED == 1
// FFT Filter Implementation using FFTW3
void rsFFTFilterFFTW(struct rsFFTFilterParams p, double *data) {
    
    /* Pad data with zeros if desired */
    double *unpaddedData;
    double *tmp;
    
    unpaddedData = data;
    data = NULL;
    data = (double*) fftw_malloc(p.paddedT*sizeof(double));
    tmp  = (double*) fftw_malloc(p.paddedT*sizeof(double));
    
    for (int i=0; i<p.T; i=i+1) {
        data[i] = unpaddedData[i];
    }
    
    for (int i=p.T; i<p.paddedT; i=i+1) {
        data[i] = 0.0;
    }
    
    /* FFT */
    fftw_execute_r2r(p.plan_r2hc, data, tmp);
    
    /* Multiply frequency bins with attenuation weight */
    for (int i = 0; i<p.paddedT; i=i+1) {
        tmp[i] = tmp[i] * p.binAttenuation[i];
    }
    
    /* Inverse FFT */
    fftw_execute_r2r(p.plan_hc2r, tmp, data);
    
    /* Remove padding and normalize */
    for (int i=0; i<p.T; i=i+1) {
        unpaddedData[i] = data[i] / p.paddedT;
    }
    
    /* Free memory */
    fftw_free(data);
    fftw_free(tmp);
    data = unpaddedData;
}
#endif

// FFT Filter Implementation using GSL
void rsFFTFilterGSL(struct rsFFTFilterParams p, double *data) {
    
    /* Prepare FFT Filtering */
    gsl_fft_real_wavetable        *real;
    gsl_fft_halfcomplex_wavetable *hc;
    gsl_fft_real_workspace        *work;
    
    work = gsl_fft_real_workspace_alloc(p.paddedT);
    real = gsl_fft_real_wavetable_alloc(p.paddedT);
    hc   = gsl_fft_halfcomplex_wavetable_alloc(p.paddedT);
    
    /* Pad data with zeros if desired */
    double *unpaddedData;
    if ( p.paddedT > p.T ) {
        unpaddedData = data;
        data = NULL;
        data = malloc(p.paddedT*sizeof(double));
        
        for (int i=0; i<p.T; i=i+1) {
            data[i] = unpaddedData[i];
        }
        
        for (int i=p.T; i<p.paddedT; i=i+1) {
            data[i] = 0.0;
        }
    }
    
    /* FFT */
    gsl_fft_real_transform(data, 1, p.paddedT, real, work);
    
    /* Multiply frequency bins with attenuation weight */
    for (int i = 0; i<p.paddedT; i=i+1) {
        data[i] = data[i] * p.binAttenuation[i];
    }
    
    /* Inverse FFT */
    gsl_fft_halfcomplex_inverse(data, 1, p.paddedT, hc, work);
    
    /* Remove padding */
    if ( p.paddedT > p.T ) {
        for (int i=0; i<p.T; i=i+1) {
            unpaddedData[i] = data[i];
        }
        free(data);
        data = unpaddedData;
    }
    
    /* Free memory */
    gsl_fft_real_wavetable_free(real);
    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);
}

void rsFFTFilter(struct rsFFTFilterParams p, double *data) {

#if RS_FFTW_ENABLED == 1
    if (rsFFTFilterEngine == RSFFTFILTER_ENGINE_FFTW) {
        rsFFTFilterFFTW(p, data);
    } else if (rsFFTFilterEngine == RSFFTFILTER_ENGINE_GSL) {
#endif
        rsFFTFilterGSL(p, data);
#if RS_FFTW_ENABLED == 1
    }
#endif
    
}

void rsFFTFilterFree(struct rsFFTFilterParams p) {
    free(p.frequencyBins);
    free(p.binAttenuation);
    
#if RS_FFTW_ENABLED == 1
    if (rsFFTFilterEngine == RSFFTFILTER_ENGINE_FFTW) {
        fftw_destroy_plan(p.plan_r2hc);
        fftw_destroy_plan(p.plan_hc2r);
    }
#endif
}

void rsFFTSetEngine(int engine) {
    
    switch (engine) {
#if RS_FFTW_ENABLED == 1
        case RSFFTFILTER_ENGINE_FFTW:
#endif
        case RSFFTFILTER_ENGINE_GSL:
            rsFFTFilterEngine = engine;
            return;
    }
    
    fprintf(stderr, "Requested FFT engine is not available!\n");
    exit(EXIT_FAILURE);
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

double rsSigmoidRolloff(const double nBins, const double rolloff, const double bin)
{
    return rsSigmoid(rolloff/nBins, bin);
}

double rsSigmoid(const double rolloff, const double x)
{
    return 1.0 - tanh(rolloff*M_PI*pow(x, 2.0));
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