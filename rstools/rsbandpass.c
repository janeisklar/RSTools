#include <stdio.h>
#include <nifti1.h>
#include <fslio.h>
#include "rsniftiutils.h"
#include "rsmathutils.h"

int show_help( void )
{
   printf(
      "rsbandpass: Given a 4D-Nifti and a frequency band this tool will apply\n"
      "            FFT filtering on it.\n"
      "\n"
   );
    
   printf(
      "basic usage:  rsbandpass -input <volume> -filtered <volume> -f1 <double> -f2 <double> -sampling_rate <double> [-savemask <mask>] [-verbose]\n"
      "\n"
   );
    
   printf(
      "options:\n"
   );

   printf(
      "   -help                   : show this help\n"
   );
 
   printf(
      "   -input <volume>         : the volume to be regressed\n"
   );
    
   printf(
      "   -filtered <volume>      : the volume in which the filtered data will be saved\n"
   );
    
   printf(
      "   -f1 <double>            : the lower frequency of the bandpass filter\n"
   );
   
   printf(
      "   -f2 <double>            : the upper frequency of the bandpass filter\n"
   );
    
   printf(
      "   -samplingrate <double> : the sampling_rate used for the FFT\n"
   );
    
   printf(
      "   -mask <mask>            : a mask specifying the ROI for improved performance\n"
   );
    
   printf(
      "   -savemask <mask>        : optional path where the rescaled mask specified with\n"
      "                             -mask will be saved. The saved file with have the same\n"
      "                             dimensions as the input volume.\n"
   );
   
   printf(
      "   -v[erbose]              : show debug information\n"
      "\n"
   );
    
   return 0;
}

void testFFTFilter();

int main(int argc, char * argv[])
{
    FSLIO *fslio;
	void *buffer;
	size_t buffsize;
	
	char *inputpath        = NULL;
	char *maskpath         = NULL;
	char *saveFilteredPath = NULL;
    char *savemaskpath     = NULL;
	
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
    
    double f1 = -1.0, f2 = -1.0, sampling_rate = -1.0;
    
    BOOL verbose = FALSE;
	
	int ac;
    
	if( argc < 2 ) return show_help();
    
    /* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if( ! strncmp(argv[ac], "-h", 2) ) {
			return show_help();
		} else if ( ! strcmp(argv[ac], "-input") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -input\n");
				return 1;
			}
			inputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-filtered") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -filtered\n");
				return 1;
			}
			saveFilteredPath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-f1") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -f1\n");
				return 1;
			}
			f1 = atof(argv[ac]);  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-f2") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -f2\n");
				return 1;
			}
			f2 = atof(argv[ac]);  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-samplingrate") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -samplingrate\n");
				return 1;
			}
			sampling_rate = atof(argv[ac]);  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-m", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m\n");
				return 1;
			}
			maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-s", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -savemask\n");
				return 1;
			}
			savemaskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-v", 2) ) {
			verbose = TRUE;
		} else if ( ! strcmp(argv[ac], "-test") ) {
            testFFTFilter();
			return 0;
		} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
		}
	}
	
	if ( inputpath == NULL ) {
		fprintf(stderr, "No input volume specified(-input)!\n");
		return 1;
	}
	
	if ( saveFilteredPath == NULL ) {
		fprintf(stderr, "An output path for the filtered data must be specified(-filtered)!\n");
		return 1;
	}
	
	if ( f1 < 0 || f2 < 0 || sampling_rate < 0 ) {
		fprintf(stderr, "Bandpass frequencies and sampling rate have to be specified!(-f1, -f2, -samplingrate)!\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Input file: %s\n", inputpath);
        fprintf(stdout, "Mask file: %s\n", maskpath);
        fprintf(stdout, "Filtered file: %s\n", saveFilteredPath);
        fprintf(stdout, "F1: %.4f\n", f1);
        fprintf(stdout, "F2: %.4f\n", f2);
        fprintf(stdout, "Sampling rate: %.4f\n", sampling_rate);
    }
    
    fslio = FslOpen(inputpath, "rb");
    if (fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n",inputpath);
        return 1;
    }
    
	/* determine dimensions */
	FslGetDim(fslio, &xDim, &yDim, &zDim, &vDim);
	   
    if ( verbose ) {
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", xDim, yDim, zDim, vDim);
    }
    
    if (fslio->niftiptr->scl_slope != 0) {
        slope = fslio->niftiptr->scl_slope;
        inter = fslio->niftiptr->scl_inter;
    }
	
	/* determine datatype and initalize buffer */
	dt = FslGetDataType(fslio, &pixtype);
	
    if ( verbose ) {
        fprintf(stdout, "Dt: %ld Pixtype: %d\n", dt, pixtype);
    }

    /* prepare filtered file */
    FSLIO *fslioFiltered;
   	void *filteredBuffer;
    fslioFiltered = FslOpen(saveFilteredPath, "wb");
        
    if (fslioFiltered == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", saveFilteredPath);
        return 1;
    }
    
    FslCloneHeader(fslioFiltered, fslio);
    FslSetDim(fslioFiltered, xDim, yDim, zDim, vDim);
    FslSetDimensionality(fslioFiltered, 4);
    FslSetDataType(fslioFiltered, pixtype);
    FslWriteHeader(fslioFiltered);
    
    buffsize = (size_t)((size_t)vDim*(size_t)dt/(size_t)8);
    filteredBuffer = malloc(buffsize);
    
    /* load mask */
    double ***mask = NULL;
    if ( maskpath != NULL ) {
        unsigned long nPoints = 0L;
        mask = d3matrix(zDim, yDim, xDim);
        Point3D *maskPoints = ReadMask(maskpath, xDim, yDim, zDim, &nPoints, savemaskpath, fslio, mask);
        if ( maskPoints == NULL) {
            fprintf(stderr, "\nError: Mask invalid.\n");
            FslClose(fslio);
            return 1;
        }
        free(maskPoints);
    }
        
    // Prepare buffer
    buffsize = (size_t)((size_t)vDim*(size_t)dt/(size_t)8);
    buffer = malloc(buffsize);
    double signal[vDim];

    /* Prepare empty timecourse */
    void *emptybuffer = malloc(buffsize);
    
    double v[vDim];
    for (short t=0; t<vDim; t=t+1) {
        v[t] = 0.0;
    }
    convertScaledDoubleToBuffer(fslioFiltered->niftiptr->datatype, emptybuffer, v, slope, inter, vDim, 1, 1, FALSE);
    
    /* Iterate over all voxels that are to be filtered */
    BOOL filterInitialized = FALSE;
    for (short z=0; z<zDim; z=z+1) {
        fprintf(stdout, "Filtering slice Z%03hd/%03hd\n", z+1, zDim);
        for (short y=0; y<yDim; y=y+1) {
            for (short x=0; x<xDim; x=x+1) {
                
                /* If it's not in the mask skip it to improve the performance */
                if (mask != NULL && mask[z][y][x] < 0.1) {
                    
                    /* set the value in the filtered data to 0 so that the nifti isn't empty */
                    if ( fslioFiltered != NULL ) {
                        rsWriteTimeSeries(fslioFiltered, emptybuffer, x, y, z, vDim);
                    }
                    
                    continue;
                }
                
                /* read out timecourse */
                FslReadTimeSeries(fslio, buffer, x, y, z, vDim);
                convertBufferToScaledDouble(signal, buffer, (long)vDim, slope, inter, fslio->niftiptr->datatype);
                
                /* apply filter */
                rsFFTFilter(signal, vDim, sampling_rate, f1, f2, filterInitialized==FALSE);
                filterInitialized = TRUE;
                
                /* write out filtered data */
                convertScaledDoubleToBuffer(fslioFiltered->niftiptr->datatype, filteredBuffer, signal, slope, inter, vDim, 1, 1, FALSE);
                rsWriteTimeSeries(fslioFiltered, filteredBuffer, x, y, z, vDim);
            }
        }
    }
    
    FslClose(fslioFiltered);
    free(fslioFiltered);
    free(filteredBuffer);
    
    if ( maskpath != NULL ) {
        free(mask);
    }
    
    free(buffer);
    FslClose(fslio);
    free(fslio);
    free(emptybuffer);
    
	return 0;
}



void testFFTFilter() {
    /* Create artifical data */
    double sampling_rate=1.8;
    int    T = 170;
    double f = 0.025;
    
    double angular_frequency = (2 * M_PI * f) * sampling_rate;
    double angular_frequency2 = (2 * M_PI * 0.09) * sampling_rate;
    double angular_frequency3 = (2 * M_PI * 0.07) * sampling_rate;
    double wave[T];
    double data[T];
    
    for (int i=0; i<T; i=i+1) {
        wave[i] = 0.0;
//        wave[i] = 1000L * sin((i+1L)*angular_frequency)  + 1000L;
//        wave[i] = 700L  * sin((i+1L)*angular_frequency2) + wave[i];
        wave[i] = 1000L  * sin((i+1L)*angular_frequency3) + wave[i];
        data[i] = wave[i];
        fprintf(stdout, "%.10f\n", data[i]);
    }
    
    rsFFTFilter(data, T, sampling_rate, 0.01, 0.04, FALSE);
    
    for (int i=0; i<T; i=i+1) {
        fprintf(stderr, "%.10f\n", data[i]);
    }
}














void extensionFFTFilter(double *data, int T, double sampling_rate) {
    
    /* extend data by mirroring the ends for a limited period */
    int T2 = 20;
    int TT = T+2*T2;
    double xdata[TT];
    for (int t=0; t<TT; t=t+1) {
        if ( t<T2 ) {
            xdata[t] = data[T2-t-1];
        } else if (t<(T+T2)) {
            xdata[t] = data[t-T2];
        } else {
            xdata[t] = data[TT-t-1];
        }
    }
    
    /* FFT Filtering */
    gsl_fft_real_wavetable        *real;
    gsl_fft_halfcomplex_wavetable *hc;
    gsl_fft_real_workspace        *work;
    
    work = gsl_fft_real_workspace_alloc(TT);
    real = gsl_fft_real_wavetable_alloc(TT);
    
    gsl_fft_real_transform(xdata, 1, TT, real, work);
    gsl_fft_real_wavetable_free(real);
    
    for (int i = 14; i<T-14; i=i+1) {
        xdata[i] = 0;
    }
    
    hc = gsl_fft_halfcomplex_wavetable_alloc(TT);
    
    gsl_fft_halfcomplex_inverse(xdata, 1, TT, hc, work);
    
    for (int i=T2; i<TT-T2; i=i+1) {
        //for (int i=0; i<TT; i=i+1) {
        //printf("%03d: %.2f\n", i+1, data[i]);
        printf("%.10f\n", xdata[i]);
    }
    
    for (int i=0; i<TT; i=i+1) {
        double F = i/(TT*sampling_rate);  //F = i*Fs/T for F = 0Hz ... Fs/2-Fs/T
        //    printf("%d: %.10f\n", i, F);
    }
    
    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);
}