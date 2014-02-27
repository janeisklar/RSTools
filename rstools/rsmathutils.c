/*******************************************************************
 *
 * rsmathutils.c
 *
 * 
 * André Hoffmann
 *******************************************************************/

#include "rsmathutils.h"

static int rsFFTFilterEngine = RSFFTFILTER_ENGINE_GSL;

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

struct rsFFTFilterParams rsFFTFilterInit(const int T, const long paddedT, const double sampling_rate, const double f1, const double f2, const int rolloff_method, const double rolloff, const BOOL keepMean, const int verbose) {
    
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
        
        if ( keepMean ) {
            attenuation[0] = 1.0;
        }
        
#if RS_FFTW_ENABLED == 1
    }
#endif
    
    if ( keepMean ) {
        attenuation[0] = 1.0;
    }
    
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
    const long double epsilon = 0.0000000001;
    
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
    const long double norm = fmaxl(sX, epsilon) * fmaxl(sY, epsilon) * Nm1;
    
    for (unsigned int i=0; i<length; i=i+1) {
        r = r + (long double)((long double)X[i]-meanX) * ((long double)Y[i]-meanY) / norm;
    }
    
    // Fisher's r-to-z transformation
    const long double half     = 0.5;
    const long double one      = 1.0;
    const long double minusone = -1.0;
    
    r=fminl(r, one-epsilon);
    r=fmaxl(r, minusone+epsilon);
    
    return half * logl( (one + r) / (one - r) );
}

double rsFastZCorrelation(const double* X, const double* Y, const size_t length)
{
    /*
     * This implementation is based on:
     * Weissenbacher, Andreas, et al.
     * "Correlations and anticorrelations in resting-state functional connectivity MRI: a quantitative comparison of preprocessing strategies."
     * Neuroimage 47.4 (2009): 1408-1416.
     */
    
    const double N = length;
    const double Nm1 = length - 1;
    const double epsilon = 0.0000000001;
    
    // compute means
    double meanX = 0.0;
    double meanY = 0.0;
    
    for (unsigned int i=0; i<length; i=i+1) {
        meanX = meanX + X[i] / N;
        meanY = meanY + Y[i] / N;
    }
    
    // compute standard scores(std. dev)
    double sX = 0.0;
    double sY = 0.0;
    
    for (unsigned int i=0; i<length; i=i+1) {
        sX = sX + pow(X[i]-meanX, 2.0) / Nm1;
        sY = sY + pow(Y[i]-meanY, 2.0) / Nm1;
    }
    
    sX = sqrt(sX);
    sY = sqrt(sY);
    
    // compute correlation coeeficient
    double r = 0.0;
    const double norm = fmax(sX, epsilon) * fmax(sY, epsilon) * Nm1;
    
    for (unsigned int i=0; i<length; i=i+1) {
        r = r + (X[i]-meanX) * (Y[i]-meanY) / norm;
    }
    
    // Fisher's r-to-z transformation
    const double half     = 0.5;
    const double one      = 1.0;
    const double minusone = -1.0;
    
    r=fmin(r, one-epsilon);
    r=fmax(r, minusone+epsilon);
    
    return half * log( (one + r) / (one - r) );
}

double rsTCorrelation(const double* X, const double* Y, const size_t length)
{
    const double r = rsCorrelation(X, Y, length);
    return r * sqrt( (length-2.0) / (1.0-r*r));
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

long double *rsFirstEigenvector(const double** A, const long n, const long maxIterations, const double precision, const BOOL verbose)
{
    // Power iteration method
    // http://en.wikipedia.org/wiki/Power_iteration
    
    long double *tmp = malloc(n*sizeof(long double));
        
    // initialize seed eigenvector
    long double *b = malloc(n*sizeof(long double));
    for (long i=0; i<n; i=i+1) {
        b[i] = 1.0;
    }
    
    BOOL converged = FALSE;
    
    // run power iteration
    for (long iteration=1L; iteration<maxIterations; iteration=iteration+1) {
                
        // calculate the matrix-by-vector product tmp=Ab
        free(tmp);
        tmp = rsMatrixByVectorProduct(A, b, n, n);
        
        // calculate the vector norm(euclidean)
        const long double norm = rsEuclideanNorm(tmp, n);
        long double invNorm = ((long double)1) / norm;
        
        // normalize it for the next iteration
        rsScaleVector(tmp, n, invNorm);
        
        // check if it has converged
        rsVectorSwap(b, tmp, n);
        rsVectorSub(tmp, b, n);
        long double mean = fabsl(rsVectorMean(tmp, n));
        
        if (verbose) {
            fprintf(stdout, "Iteration %ld/%ld Mean Difference: %.10Lf Target: %.10f\n", iteration, maxIterations, mean, precision);
        }
        
        if ( mean <= precision ) {
            break;
        }
    }
    
    free(tmp);
    
    return b;
}

// Originally taken from: https://gist.github.com/microo8/4065693
struct rsPCAResult rsGenericPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose)
{
    /*
    @param data - matrix of data vectors, MxN matrix, each column is a data vector, M - dimension, N - data vector count
    @param minVariance - min percentage of variance that should be retained
    @param nComponents - number of components that will be returned. ignored if less than 1
    */
	
    assert(data != NULL);
    unsigned int i;
    unsigned int j;
    unsigned int rows = data->size1;
    unsigned int cols = data->size2;
	struct rsPCAResult result;

	if ( verbose ) {
		fprintf(stdout, "Running PCA on a %dx%d matrix\n", rows, cols);
	}

    gsl_vector* mean = gsl_vector_alloc(rows);
 
    for(i = 0; i < rows; i++) {
        gsl_vector_set(mean, i, gsl_stats_mean(data->data + i * cols, 1, cols));
    }
 
    // Get mean-substracted data into matrix mean_substracted_data.
    gsl_matrix* mean_substracted_data = gsl_matrix_alloc(rows, cols);
    gsl_matrix_memcpy(mean_substracted_data, data);
    for(i = 0; i < cols; i++) {
        gsl_vector_view mean_substracted_point_view = gsl_matrix_column(mean_substracted_data, i);
        gsl_vector_sub(&mean_substracted_point_view.vector, mean);
    }
    gsl_vector_free(mean);
 
    // Compute Covariance matrix
    gsl_matrix* covariance_matrix = gsl_matrix_alloc(rows, rows);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(cols - 1), mean_substracted_data, mean_substracted_data, 0.0, covariance_matrix);
    gsl_matrix_free(mean_substracted_data);
 
    // Get eigenvectors, sort by eigenvalue.
    gsl_vector* eigenvalues = gsl_vector_alloc(rows);
    gsl_matrix* eigenvectors = gsl_matrix_alloc(rows, rows);
    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(rows);
    gsl_eigen_symmv(covariance_matrix, eigenvalues, eigenvectors, workspace);
    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(covariance_matrix);
 
    // Sort the eigenvectors
    gsl_eigen_symmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);

	if ( verbose ) {
		fprintf(stdout, "\nPCAs(eigenvectors):\n");
		
		for ( i = 0; i<rows; i=i+1 ) {
			for ( j = 0; j<rows; j=j+1) {
            	fprintf(stdout, "% 5.10f\t", gsl_matrix_get(eigenvectors, i, j));
			}
			fprintf(stdout, "\n");
        }
	}

	// Determine how many components to keep
	int L = 0;
	double explainedVariance = 0.0;
	double totalVariance = 0.0;
	for (unsigned int i=0; i<rows; i=i+1)
		totalVariance = totalVariance + gsl_vector_get(eigenvalues, i);
	
	if ( verbose ) {
		fprintf(stdout, "\nTotal variance: %.10f\n", totalVariance);
		fprintf(stdout, "\nEigenvalues:\n");
		for ( j = 0; j<rows; j=j+1) {
        	fprintf(stdout, "%.10f\t", gsl_vector_get(eigenvalues, j)/totalVariance);
		}
		fprintf(stdout, "\n\n");
		fprintf(stdout, "Determining the number of components to keep:\n");
	}
		
	do {
		explainedVariance = explainedVariance + gsl_vector_get(eigenvalues, L);
		
		if ( verbose )
			fprintf(stdout, "% 3d %.10f\n", L, (explainedVariance/totalVariance));
			
		L = L + 1;
	} while( L < nComponents || ((explainedVariance/totalVariance) < minVariance && nComponents < 1) );
	
	if ( verbose )
		fprintf(stdout, "\n");
	
 	
    gsl_matrix_view L_eigenvectors = gsl_matrix_submatrix(eigenvectors, 0, 0, rows, L);
    gsl_vector_view L_eigenvalues  = gsl_vector_subvector(eigenvalues, 0, L);

    // Project the original dataset
    result.transformed = gsl_matrix_alloc(L, cols); // transformed is a n LxN matrix, each column is the original data vector with reduced dimension from M to L
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &L_eigenvectors.matrix, data, 0.0, result.transformed);

	// Extract the first L eigenvectors
 	result.eigenvectors = gsl_matrix_alloc(L, rows);
	gsl_matrix_transpose_memcpy(result.eigenvectors, &(L_eigenvectors.matrix));
	
	// Extract the first L eigenvalues
	result.eigenvalues = gsl_vector_alloc(L);
	gsl_vector_memcpy(result.eigenvalues, &(L_eigenvalues.vector));

	// Also store all eigenvalues
	result.eigenvalues_all = eigenvalues;
	
	gsl_matrix_free(eigenvectors);
	
	return result;
}

struct rsPCAResult rsPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose)
{
	return rsGenericPCA(data, minVariance, nComponents, verbose);
}

struct rsPCAResult rsTPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose)
{
	gsl_matrix* dataT = gsl_matrix_alloc(data->size2, data->size1);
	gsl_matrix_transpose_memcpy(dataT, data);
	
	struct rsPCAResult result = rsGenericPCA(dataT, minVariance, nComponents, verbose);
	gsl_matrix_free(dataT);
	
	return result;
}

void rsPCAResultFree(struct rsPCAResult result)
{
	gsl_matrix_free(result.transformed);
	gsl_matrix_free(result.eigenvectors);
	gsl_vector_free(result.eigenvalues);
	gsl_vector_free(result.eigenvalues_all);
}

struct rsCTPResult rsCTP(const gsl_matrix* A, const gsl_matrix* B, int nComponents, BOOL verbose)
{
	
    assert(A != NULL);
    assert(B != NULL);
	assert(A->size1 == B->size1);

    long i, j;
    unsigned int rows  = A->size1;
    unsigned int colsA = A->size2;
    unsigned int colsB = B->size2;
	struct rsCTPResult result;

	if ( verbose ) {
		fprintf(stdout, "Running Common Temporal Patterns Analysis on a %dx%d and %dx%d matrix\n", rows, colsA, rows, colsB);
	}
	
	// compute means
    gsl_vector* meanA = gsl_vector_alloc(rows);
    gsl_vector* meanB = gsl_vector_alloc(rows);
 	
	#pragma omp parallel num_threads(rsGetThreadsNum()) shared(colsA,colsB,rows) private(i,j)
	{
        #pragma omp for schedule(guided)
       	for(i = 0; i < rows; i=i+1) {
			double mean = 0.0;
			for(j=0; j<colsA; j=j+1) {
				mean = mean + gsl_matrix_get(A,i,j) / colsA;
			}
			gsl_vector_set(meanA, i, mean);
			
			mean = 0.0;
			for(j=0; j<colsB; j=j+1) {
				mean = mean + gsl_matrix_get(B,i,j) / colsB;
			}
        	gsl_vector_set(meanB, i, mean);
		}
    }
 
    // Get mean-substracted data into matrix demeanedA and demeanedB.
    gsl_matrix* demeanedA;
    gsl_matrix* demeanedB;
    
	#pragma omp parallel sections num_threads(rsGetThreadsNum())
	{
		#pragma omp section
		{
			demeanedA = gsl_matrix_alloc(rows, colsA);
			gsl_matrix_memcpy(demeanedA, A);
		}
		
		#pragma omp section
		{
			demeanedB = gsl_matrix_alloc(rows, colsB);
    		gsl_matrix_memcpy(demeanedB, B);
		}
	}

	#pragma omp parallel num_threads(rsGetThreadsNum()) shared(colsA,colsB) private(i)
	{
        #pragma omp for schedule(guided)
    	for(i = 0; i < colsA; i++) {
        	gsl_vector_view demeanedPointViewA = gsl_matrix_column(demeanedA, i);
        	gsl_vector_sub(&demeanedPointViewA.vector, meanA);
		}
		
        #pragma omp for schedule(guided)
    	for(i = 0; i < colsB; i++) {
        	gsl_vector_view demeanedPointViewB = gsl_matrix_column(demeanedB, i);
        	gsl_vector_sub(&demeanedPointViewB.vector, meanB);
		}
    }
 
    // Compute Covariance matrices
	gsl_matrix* covA;
	gsl_matrix* covB;
	
	#pragma omp parallel sections num_threads(rsGetThreadsNum())
	{
		#pragma omp section
		{
		    gsl_vector_free(meanA);
			covA = gsl_matrix_alloc(rows, rows);
			gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(colsA - 1), demeanedA, demeanedA, 0.0, covA);
    		gsl_matrix_free(demeanedA);
		}
	
		#pragma omp section
		{
		    gsl_vector_free(meanB);
			covB = gsl_matrix_alloc(rows, rows);
    		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(colsB - 1), demeanedB, demeanedB, 0.0, covB);
    		gsl_matrix_free(demeanedB);
		}
	}
	
	gsl_matrix_add(covB, covA);
 
    // Get eigenvectors, sort by eigenvalue.
    gsl_vector* eigenvalues = gsl_vector_alloc(rows);
    gsl_matrix* eigenvectors = gsl_matrix_alloc(rows, rows);
    gsl_eigen_gensymmv_workspace* workspace = gsl_eigen_gensymmv_alloc(rows);
    gsl_eigen_gensymmv(covA, covB, eigenvalues, eigenvectors, workspace);
    gsl_eigen_gensymmv_free(workspace);
    gsl_matrix_free(covA);
    gsl_matrix_free(covB);
 
    // Sort the eigenvectors
    gsl_eigen_gensymmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);

	// Extract the first and last nComponents eigenvectors
	int L = 2 * nComponents;
	
    gsl_matrix* L_eigenvectors = gsl_matrix_alloc(rows, L);
    gsl_vector* L_eigenvalues  = gsl_vector_alloc(L);
	
	#pragma omp parallel num_threads(rsGetThreadsNum()) shared(L,rows) private(i)
	{
        #pragma omp for schedule(guided)
		for ( i = 0; i<nComponents; i=i+1 ) {
			const unsigned i2 = L-i-1;
			gsl_vector_view viewCol  = gsl_matrix_column(eigenvectors, i);
			gsl_vector_view viewCol2 = gsl_matrix_column(eigenvectors, rows-i-1);
			gsl_matrix_set_col(L_eigenvectors, i,  &(viewCol.vector));
			gsl_matrix_set_col(L_eigenvectors, i2, &(viewCol2.vector));
			gsl_vector_set(L_eigenvalues, i, gsl_vector_get(eigenvalues, i));
			gsl_vector_set(L_eigenvalues, i2, gsl_vector_get(eigenvalues, rows-i-1));
		}
	}
 	
	#pragma omp parallel sections num_threads(rsGetThreadsNum())
	{
		// Project the original dataset
		#pragma omp section
		{
			result.transformedA = gsl_matrix_alloc(L, colsA);
		    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L_eigenvectors, A, 0.0, result.transformedA);
		}
		
		#pragma omp section
		{
    		result.transformedB = gsl_matrix_alloc(L, colsB);
			gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L_eigenvectors, B, 0.0, result.transformedB);
		}
		
		// Store eigenvectors
		#pragma omp section
		{
			result.eigenvectors = gsl_matrix_alloc(L, rows);
			gsl_matrix_transpose_memcpy(result.eigenvectors, L_eigenvectors);
		}
		
		// Store eigenvalues
		#pragma omp section
		{
			result.eigenvalues     = L_eigenvalues;
			result.eigenvalues_all = eigenvalues;
		}
	}

	
	return result;
}

void rsCTPResultFree(struct rsCTPResult result)
{
	gsl_matrix_free(result.transformedA);
	gsl_matrix_free(result.transformedB);
	gsl_matrix_free(result.eigenvectors);
	gsl_vector_free(result.eigenvalues);
	gsl_vector_free(result.eigenvalues_all);
}

double **d2matrix(int yh, int xh)
{
    long int j;
    long int nrow = yh+1;
    long int ncol = xh+1;
    double **t;
    
    
    /** allocate pointers to ydim */
    t=(double **) malloc((size_t)((nrow)*sizeof(double*)));
    if (!t) RSIOERR("d2matrix: allocation failure");
    
    /** allocate pointers for zdim */
    t[0]=(double *) malloc((size_t)((ncol*nrow)*sizeof(double)));
    if (!t[0]) RSIOERR("d2matrix: allocation failure");
    
    /** point everything to the data blob */
    for(j=1L;j<nrow;j++) t[j]=t[j-1]+ncol;
    
    return t;
}

/*
 * Computes the matrix by vector product
 * y = A*x
 * n..width of the matrix and vector length
 * m..height of the matrix
 * The matrix is assumed to be in the following format A[m][n]
 */
long double *rsMatrixByVectorProduct(const double **A, const long double *x, const long n, const long m) {
    long double *y = malloc(m*sizeof(long double));
    long i,j;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,j) shared(y)
    {
        #pragma omp for schedule(guided)
        for (i=0; i<m; i=i+1L) {
            long double product = 0.0;
            
            for (j=0; j<n; j=j+1) {
                product += (long double)A[i][j] * (long double)x[j];
            }
            
            y[i] = (double)product;
        }
    }
    
    return y;
}

long double rsEuclideanNorm(const long double *x, const long n)
{
    long double tss = 0.0;
    
    for (long i=0; i<n; i=i+1) {
        tss += x[i];
    }
    
    return sqrtl(tss);
}

void rsScaleVector(long double *x, const long n, const long double factor)
{
    long i;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i) shared(x)
    {
        #pragma omp for schedule(guided)
        for (i=0L; i<n; i=i+1L) {
            x[i] *= factor;
        }
    }
}

void rsVectorSub(long double *x, const long double *y, const long n)
{
    long i;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i) shared(x,y)
    {
        #pragma omp for schedule(guided)
        for (i=0L; i<n; i=i+1L) {
            x[i] -= y[i];
        }
    }
}

void rsVectorSwap(long double *x, long double *y, const long n)
{
    long double tmp;
    long i;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,tmp) shared(x,y)
    {
        #pragma omp for schedule(guided)
        for (i=0L; i<n; i=i+1L) {
            tmp  = x[i];
            x[i] = y[i];
            y[i] = tmp;
        }
    }
}

BOOL rsVectorContains(const long *x, const long n, const long element) {
	int i;
	for ( i=0; i<n; i=i+1 ) {
		if ( x[i] == element ) {
			return TRUE;
		}
	}
	
	return FALSE;
}

void rsMatrixConversion(double **A, const long m, const long n, const int mode)
{
    long i,j;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,j) shared(A)
    {
        #pragma omp for schedule(guided)
        for (i=0; i<m; i=i+1L) {
            if ( mode == RSMATRIXCONVERSION_ABSOLUTE ) {
                for (j=0; j<n; j=j+1) {
                    A[i][j] = fabs(A[i][j]);
                }
            } else if ( mode == RSMATRIXCONVERSION_POSITIVE ) {
                for (j=0; j<n; j=j+1) {
                    if ( A[i][j] < 0.0 ) {
                        A[i][j] = 0.0;
                    }
                }
            } else if ( mode == RSMATRIXCONVERSION_NEGATIVE ) {
                for (j=0; j<n; j=j+1) {
                    if ( A[i][j] >= 0.0 ) {
                        A[i][j] = 0.0;
                    } else {
                        A[i][j] = fabs(A[i][j]);
                    }
                }
            } else if ( mode == RSMATRIXCONVERSION_SCALED ) {
                // This is used to shift z-transformed correlations so that they're all positive.
                // Therefore the smallest possible(as implemented) correlation coefficient is computed.
                const double epsilon = 0.0000000001;
                const double shiftingConstant = fabs(log(epsilon / (2.0 - epsilon))) + epsilon;
                
                for (j=0; j<n; j=j+1) {
                    A[i][j] = A[i][j] + shiftingConstant;
                }
            }
        }
    }
}

long double rsVectorMean(const long double *x, const long n)
{
    long double mean = 0.0;
    
    for (long i=0L; i<n; i=i+1L) {
        mean += x[i];
    }
    
    return (double) (mean/(unsigned long)n);
}

void rs_matrix_fprintf(FILE *stream, const double **A, const long m, const long n, const char* fmt)
{
    for (long i=0; i<m; i=i+1) {
        for (long j=0; j<n; j=j+1) {
            fprintf(stream, fmt, A[i][j]);
            fprintf(stream, " ");
        }
        fprintf(stream, "\n");
    }
}

void rs_vector_fprintf(FILE *stream, const double *x, const long n, const char* fmt)
{
    for (long j=0; j<n; j=j+1) {
        fprintf(stream, fmt, x[j]);
        fprintf(stream, " ");
    }
    fprintf(stream, "\n");
}

void rs_vector_fprintfl(FILE *stream, const long double *x, const long n, const char* fmt)
{
        for (long j=0; j<n; j=j+1) {
            fprintf(stream, fmt, x[j]);
            fprintf(stream, " ");
        }
        fprintf(stream, "\n");
}

int rs_gsl_matrix_fprintf(FILE *stream,gsl_matrix *m,char *fmt)
{
    size_t rows=m->size1;
    size_t cols=m->size2;
    size_t row,col,ml;
    int fill;
    char buf[100];
    gsl_vector *maxlen;
    
    maxlen=gsl_vector_alloc(cols);
    for (col=0;col<cols;++col) {
        ml=0;
        for (row=0;row<rows;++row) {
            sprintf(buf,fmt,gsl_matrix_get(m,row,col));
            if (strlen(buf)>ml)
                ml=strlen(buf);
        }
        gsl_vector_set(maxlen,col,ml);
    }
    
    for (row=0;row<rows;++row) {
        for (col=0;col<cols;++col) {
            sprintf(buf,fmt,gsl_matrix_get(m,row,col));
            fprintf(stream,"%s",buf);
            fill=gsl_vector_get(maxlen,col)+1-strlen(buf);
            while (--fill>=0)
                fprintf(stream," ");
        }
        fprintf(stream,"\n");
    }
    gsl_vector_free(maxlen);
    return 0;
}

BOOL rsSaveMatrix(const char *filename, const double** A, const long m, const long n)
{
    FILE *file;
    file = fopen(filename, "wb");
    BOOL success = TRUE;
    const size_t blockSize = sizeof(double);
    
    for (long row=0L; row<m && success; row=row+1) {
        success = success && fwrite(A[row], blockSize, n, file) == n;
    }
    fclose(file);
    
    return success;
}

BOOL rsLoadMatrix(const char *filename, double** A, const long m, const long n)
{
    FILE *file;
    file = fopen(filename, "r");
    BOOL success = TRUE;
    const size_t blockSize = sizeof(double);
    
    for (long row=0L; row<m && success; row=row+1) {
        success = success && fread(A[row], blockSize, n, file) == n;
    }
    fclose(file);
    
    return success;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// See: http://en.wikipedia.org/wiki/Student's_t-test#One-sample_t-test                                     //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double rsOneSampleTTest(const double *data, const unsigned int length, const double mu)
{
    // compute mean
    double mean = 0.0;
    
    for (unsigned int i=0; i<length; i=i+1) {
        mean = mean + data[i] / length;
    }
    
    // compute standard scores(std. dev)
    double std = 0.0;
    
    for (unsigned int i=0; i<length; i=i+1) {
        std = std + pow(data[i]-mean, 2.0) / ((double)length-1.0);
    }
    
    std = sqrt(std);
    
    return (mean - mu) / (std / sqrt(length));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following two functions are taken from:                                                              //
// http://www.spraak.org/documentation/doxygen/src/lib/math/erfinv.c/view                                   //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*************************************************************
 *  __  _   _   _   _                                        *
 * (_  |_) |_) |_| |_| |_/  Speech Processing, Recognition & *
 * __) |   | \ | | | | | \  Automatic Annotation Kit         *
 *                                                           *
 * Copyright 2006,2007,2008 K.U.Leuven                       *
 *                                                           *
 * Use of this software is governed by a License Agreement.  *
 * Modifications should be properly credited in line with    *
 * the License Agreement.                                    *
 *************************************************************/

/*
 * Calculate the inverse error function.
 * The fast version is only correct up to 6 digits.
 */
float rsFastErfInv(float x) {
    float tmp;
    int   neg;
    
    if((neg=(x < 0.0)))
        x = -x;
    if(x <= 0.7) {
        tmp = x*x;
        x *= (((-0.140543331*tmp+0.914624893)*tmp-1.645349621)*tmp+0.886226899)/((((0.012229801*tmp-0.329097515)*tmp+1.442710462)*tmp-2.118377725)*tmp+1.0);
    } else {
        tmp = sqrt(-log(0.5*(1.0-x)));
        x = (((1.641345311*tmp+3.429567803)*tmp-1.624906493)*tmp-1.970840454)/((1.637067800*tmp+3.543889200)*tmp+1.0);
    }
    return(neg?-x:x);
}

/*
 * Calculate the inverse error function.
 * Uses fast_erfinv and performs a two steps of Newton-Raphson correction to
 * achieve full accuracy.
 */
double rsErfInv(const double x) {
    const double tmp = rsFastErfInv(x);
    return tmp - (erf(tmp)-x)*exp(tmp*tmp)*0.886226925452757941 - (erf(tmp)-x)*exp(tmp*tmp)*0.886226925452757941;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following two functions are taken from:                                                              //
// https://gist.github.com/timflutre/1784199                                                                //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Resolve a sequence of ties.
 *
 * The input ranks array is expected to take the same value for all indices in
 * tiesTrace. The common value is recoded with the average of the indices. For
 * example, if ranks = <5,8,2,6,2,7,1,2> and tiesTrace = <2,4,7>, the result
 * will be <5,8,3,6,3,7,1,3>.
 *
 * Source: http://commons.apache.org/math/apidocs/src-html/org/apache/commons/math/stat/ranking/NaturalRanking.html#line.312
 */
void rsRankingResolveTies(double *ranks, const size_t *tiesTrace, const size_t n_ties) {
    size_t i;
    
    // constant value of ranks over tiesTrace
    const double c = ranks[tiesTrace[0]];
    
    // new rank (ie. the average of the current indices)
    const double avg = (2*c + n_ties - 1) / 2;
    
    for(i=0; i<n_ties; ++i) {
        ranks[tiesTrace[i]] = avg;
    }
}

/*
 * Rank data using the natural ordering on doubles, ties being resolved by taking their average.
 *
 * Source: http://commons.apache.org/math/apidocs/src-html/org/apache/commons/math/stat/ranking/NaturalRanking.html#line.190
 */
void rsSpearmannRank(double *ranks, const double *data, const size_t n) {
    size_t i;
    double *d = malloc(sizeof(double)*n);
    size_t *p = malloc(sizeof(size_t)*n);
    
    // copy the input data and sort them
    for(i=0; i<n; ++i)
        d[i] = data[i];
    gsl_sort(d, 1, n);
    
    // get the index of the input data as if they were sorted
    gsl_sort_index(p, data, 1, n);
    
    // walk the sorted array, filling output array using sorted positions, resolving ties as we go
    size_t pos = 1;
    ranks[p[0]] = pos;
    size_t n_ties = 1;
    size_t * tiesTrace = (size_t*) calloc (1, sizeof(size_t));
    tiesTrace[0] = p[0];
    for(i=1; i<n; ++i) {
        if(d[i] - d[i-1] > 0) {
            pos = i + 1;
            if(n_ties > 1) {
                rsRankingResolveTies(ranks, tiesTrace, n_ties);
            }
            tiesTrace = (size_t*) realloc(tiesTrace, sizeof(size_t));
            n_ties = 1;
            tiesTrace[0] = p[i];
        } else {
            ++n_ties;
            tiesTrace = (size_t*) realloc(tiesTrace, n_ties * sizeof(size_t));
            tiesTrace[n_ties-1] = p[i];
        }
        ranks[p[i]] = pos;
    }
    if(n_ties > 1) {
        rsRankingResolveTies(ranks, tiesTrace, n_ties);
    }
    
    free(tiesTrace);
    free(d);
    free(p);
}