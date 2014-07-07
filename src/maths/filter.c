#include "filter.h"

static int rsFFTFilterEngine = RSFFTFILTER_ENGINE_GSL;

rsFFTFilterParams* rsFFTFilterInit(const int T, const long paddedT, const double sampling_rate, const double f1, const double f2, const int rolloff_method, const double rolloff, const BOOL keepMean, const int verbose) {
    
    rsFFTFilterParams *p = rsMalloc(sizeof(rsFFTFilterParams));
    double *F = rsMalloc(paddedT * sizeof(double));
    double *attenuation = rsMalloc(paddedT*sizeof(double));

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
        
        p->plan_r2hc = fftw_plan_r2r_1d(paddedT, in, out, FFTW_R2HC, FFTW_PATIENT);
        p->plan_hc2r = fftw_plan_r2r_1d(paddedT, out, in, FFTW_HC2R, FFTW_PATIENT);

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
    
    p->frequencyBins  = F;
    p->binAttenuation = attenuation;
    p->f1             = f1;
    p->f2             = f2;
    p->verbose        = verbose;
    p->sampling_rate  = sampling_rate;
    p->T              = T;
    p->paddedT        = paddedT;
    p->rolloff_method = rolloff_method;
    p->rolloff        = rolloff;
    
    return p;
}

#if RS_FFTW_ENABLED == 1
// FFT Filter Implementation using FFTW3
void rsFFTFilterFFTW(rsFFTFilterParams *p, double *data) {
    
    /* Pad data with zeros if desired */
    double *unpaddedData;
    double *tmp;
    
    unpaddedData = data;
    data = NULL;
    data = (double*) fftw_malloc(p->paddedT*sizeof(double));
    tmp  = (double*) fftw_malloc(p->paddedT*sizeof(double));
    
    for (int i=0; i<p->T; i=i+1) {
        data[i] = unpaddedData[i];
    }
    
    for (int i=p->T; i<p->paddedT; i=i+1) {
        data[i] = 0.0;
    }
    
    /* FFT */
    fftw_execute_r2r(p->plan_r2hc, data, tmp);
    
    /* Multiply frequency bins with attenuation weight */
    for (int i = 0; i<p->paddedT; i=i+1) {
        tmp[i] = tmp[i] * p->binAttenuation[i];
    }
    
    /* Inverse FFT */
    fftw_execute_r2r(p->plan_hc2r, tmp, data);
    
    /* Remove padding and normalize */
    for (int i=0; i<p->T; i=i+1) {
        unpaddedData[i] = data[i] / p->paddedT;
    }
    
    /* Free memory */
    fftw_free(data);
    fftw_free(tmp);
    data = unpaddedData;
}
#endif

// FFT Filter Implementation using GSL
void rsFFTFilterGSL(rsFFTFilterParams *p, double *data) {
    
    /* Prepare FFT Filtering */
    gsl_fft_real_wavetable        *real;
    gsl_fft_halfcomplex_wavetable *hc;
    gsl_fft_real_workspace        *work;
    
    work = gsl_fft_real_workspace_alloc(p->paddedT);
    real = gsl_fft_real_wavetable_alloc(p->paddedT);
    hc   = gsl_fft_halfcomplex_wavetable_alloc(p->paddedT);
    
    /* Pad data with zeros if desired */
    double *unpaddedData;
    if ( p->paddedT > p->T ) {
        unpaddedData = data;
        data = NULL;
        data = rsMalloc(p->paddedT*sizeof(double));
        
        for (int i=0; i<p->T; i=i+1) {
            data[i] = unpaddedData[i];
        }
        
        for (int i=p->T; i<p->paddedT; i=i+1) {
            data[i] = 0.0;
        }
    }
    
    /* FFT */
    gsl_fft_real_transform(data, 1, p->paddedT, real, work);
    
    /* Multiply frequency bins with attenuation weight */
    for (int i = 0; i<p->paddedT; i=i+1) {
        data[i] = data[i] * p->binAttenuation[i];
    }
    
    /* Inverse FFT */
    gsl_fft_halfcomplex_inverse(data, 1, p->paddedT, hc, work);
    
    /* Remove padding */
    if ( p->paddedT > p->T ) {
        for (int i=0; i<p->T; i=i+1) {
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

void rsFFTFilter(rsFFTFilterParams *p, double *data) {

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

void rsFFTFilterFree(rsFFTFilterParams *p) {
	if ( p == NULL ) {
		return;
	}
	
    free(p->frequencyBins);
    free(p->binAttenuation);
    
#if RS_FFTW_ENABLED == 1
    if (rsFFTFilterEngine == RSFFTFILTER_ENGINE_FFTW) {
        fftw_destroy_plan(p->plan_r2hc);
        fftw_destroy_plan(p->plan_hc2r);
    }
#endif
	free(p);
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
