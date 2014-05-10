//
//  rsbandpass_common.c
//  rstools
//
//  Created by André Hoffmann on 6/28/13.
//
//

#include "rsbandpass_common.h"


void rsBandpassPrintHelp() {

    printf(
      "basic usage:  rsbandpass -input <volume> -filtered <volume> -f1 <double> -f2 <double> -sampling_rate <double> [-savemask <mask>] [-verbose]\n"
      "\n"
    );

    printf(
      "options:\n"
    );

    printf(
      "   -input <volume>          : the volume to be regressed\n"
    );

    printf(
      "   -filtered <volume>       : the volume in which the filtered data will be saved\n"
    );

    printf(
      "   -f1 <double>             : the lower frequency of the bandpass filter\n"
    );

    printf(
      "   -f2 <double>             : the upper frequency of the bandpass filter\n"
    );

    printf(
      "   -samplingrate <double>   : the sampling_rate used for the FFT\n"
    );

    printf(
      "   -mask <mask>             : a mask specifying the ROI for improved performance\n"
    );

    printf(
      "   -savemask <mask>         : optional path where the rescaled mask specified with\n"
      "                              -mask will be saved. The saved file with have the same\n"
      "                              dimensions as the input volume.\n"
    );
    
    printf(
      "   -keepMean                : keeps the first bin of the FFT(the mean) independent of\n"
      "                              tbe selected frequency range\n"
    );

    printf(
      "   -threads <int>           : number of threads used for processing\n"
    );
    
#if RS_FFTW_ENABLED == 1
    printf(
      "   -fftw                    : use FFTW3 instead of GSL for FFT\n"
    );
#endif
    
    printf(
      "   -saveattenuation <file>  : save txt file that contains the bin's frequencies and the \n"
      "                              corresponding attenuation weight that was used.\n"
    );
    
    printf(
      "   -sigmoidrolloff <double> : uses a sigmoid function for rolling of the passband. The\n"
      "                              specified number controls how fast it is rolled of with\n"
      "                              higher numbers corresponding to a quicker rolloff. A good\n"
      "                              starting point would be 10, then double-check by saving\n"
      "                              the attenuation file.\n"
    );

    printf(
      "   -v[erbose]               : show debug information\n"
      "\n"
    );
}

struct rsBandpassParameters rsBandpassInitParameters() {
    struct rsBandpassParameters p;
    
    p.inputpath            = NULL;
    p.maskpath             = NULL;
    p.savemaskpath         = NULL;
    p.saveFilteredPath     = NULL;
    p.saveAttenuationPath  = NULL;
    p.xDim                 = 0;
    p.yDim                 = 0;
    p.zDim                 = 0;
    p.vDim                 = 0;
    p.paddedT              = 0;
    p.pixtype              = 0;
    p.dt                   = 4;
    p.inter                = 0.0;
    p.slope                = 1.0;
    p.f1                   = -1.0;
    p.f2                   = -1.0;
    p.sampling_rate        = -1.0;
    p.verbose              = FALSE;
    p.fslio                = NULL;
    p.parametersValid      = FALSE;
    p.mask                 = NULL;
    p.threads              = 1;
    p.rolloff_method       = RSFFTFILTER_CUTOFF;
    p.rolloff              = 10.0;
    p.keepMean             = FALSE;
    
    return p;
}

struct rsBandpassParameters rsBandpassLoadParams(int argc, char * argv[]) {

    struct rsBandpassParameters p = rsBandpassInitParameters();
    
    /* parse parameters */
	for( int ac = 1; ac < argc; ac++ ) {
		if ( ! strcmp(argv[ac], "-input") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -input\n");
				return p;
			}
			p.inputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-filtered") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -filtered\n");
				return p;
			}
			p.saveFilteredPath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-f1") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -f1\n");
				return p;
			}
			p.f1 = atof(argv[ac]);  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-f2") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -f2\n");
				return p;
			}
			p.f2 = atof(argv[ac]);  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-samplingrate") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -samplingrate\n");
				return p;
			}
			p.sampling_rate = atof(argv[ac]);  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-m", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m\n");
				return p;
			}
			p.maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-savemask") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -savemask\n");
				return p;
			}
			p.savemaskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-threads") ) {
    		if( ++ac >= argc ) {
    			fprintf(stderr, "** missing argument for -threads\n");
    			return p;
    		}
    		p.threads = atoi(argv[ac]);
    	} else if ( ! strcmp(argv[ac], "-datalength") ) {
    		if( ++ac >= argc ) {
    			fprintf(stderr, "** missing argument for -datalength\n");
    			return p;
    		}
    		p.paddedT = atol(argv[ac]);
    	} else if ( ! strcmp(argv[ac], "-saveattenuation") ) {
            if( ++ac >= argc ) {
    			fprintf(stderr, "** missing argument for -saveattenuation\n");
    			return p;
    		}
    		p.saveAttenuationPath = argv[ac];
        } else if ( ! strcmp(argv[ac], "-sigmoidrolloff") ) {
            if( ++ac >= argc ) {
    			fprintf(stderr, "** missing argument for -sigmoidrolloff\n");
    			return p;
    		}
    		p.rolloff_method = RSFFTFILTER_SIGMOID;
            p.rolloff = atof(argv[ac]);
        } else if ( ! strcmp(argv[ac], "-fftw") ) {
            rsFFTSetEngine(RSFFTFILTER_ENGINE_FFTW);
        } else if ( ! strncmp(argv[ac], "-v", 2) ) {
			p.verbose = TRUE;
		} else if ( ! strcmp(argv[ac], "-keepMean") ) {
			p.keepMean = TRUE;
		} else if ( ! strcmp(argv[ac], "-test") ) {
            testFFTFilter();
			return p;
		} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
		}
	}
	
	if ( p.inputpath == NULL ) {
		fprintf(stderr, "No input volume specified(-input)!\n");
		return p;
	}
	
	if ( p.saveFilteredPath == NULL ) {
		fprintf(stderr, "An output path for the filtered data must be specified(-filtered)!\n");
		return p;
	}
	
	if ( p.f1 < 0 || p.f2 < 0 || p.sampling_rate < 0 ) {
		fprintf(stderr, "Bandpass frequencies and sampling rate have to be specified!(-f1, -f2, -samplingrate)!\n");
		return p;
	}
    
    if ( p.verbose ) {
        fprintf(stdout, "Input file: %s\n", p.inputpath);
        fprintf(stdout, "Mask file: %s\n", p.maskpath);
        fprintf(stdout, "Filtered file: %s\n", p.saveFilteredPath);
        fprintf(stdout, "F1: %.4f\n", p.f1);
        fprintf(stdout, "F2: %.4f\n", p.f2);
        fprintf(stdout, "Sampling rate: %.4f\n", p.sampling_rate);
    }
    
    p.fslio = FslOpen(p.inputpath, "rb");
    if (p.fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n",p.inputpath);
        return p;
    }
    
	/* determine dimensions */
	FslGetDim(p.fslio, &p.xDim, &p.yDim, &p.zDim, &p.vDim);
	   
    if ( p.verbose ) {
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", p.xDim, p.yDim, p.zDim, p.vDim);
    }

    if ( p.paddedT == 0 ) {
        p.paddedT = p.vDim;
    } else if ( p.verbose ) {
        fprintf(stdout, "Padding data to have a sampling length of %ld\n", p.paddedT);
    }
    
    if ( p.vDim > p.paddedT ) {
        fprintf(stderr, "\nError, datalength(%ld) needs to be longer or equal to the temporal length(%d) of the supplied nifti file: %s.\n",p.paddedT, p.vDim, p.inputpath);
        return p;
    }
    
    if (p.fslio->niftiptr->scl_slope != 0) {
        p.slope = p.fslio->niftiptr->scl_slope;
        p.inter = p.fslio->niftiptr->scl_inter;
    }
	
	/* determine datatype and initalize buffer */
	p.dt = FslGetDataType(p.fslio, &p.pixtype);
	
    if ( p.verbose ) {
        fprintf(stdout, "Dt: %ld Pixtype: %d\n", p.dt, p.pixtype);
    }
    
    if ( p.maskpath != NULL ) {
        unsigned long nPoints = 0L;
        p.mask = d3matrix(p.zDim, p.yDim, p.xDim);
        Point3D *maskPoints = ReadMask(p.maskpath, p.xDim, p.yDim, p.zDim, &nPoints, p.savemaskpath, p.fslio, p.mask);
        if ( maskPoints == NULL) {
            fprintf(stderr, "\nError: Mask invalid.\n");
            FslClose(p.fslio);
            return p;
        }
        free(maskPoints);
    }
    
    // Prepare FFT filter
    p.fftParams = rsFFTFilterInit(p.vDim, p.paddedT, p.sampling_rate, p.f1, p.f2, p.rolloff_method, p.rolloff, p.keepMean, p.verbose);
    
    if ( p.saveAttenuationPath != NULL ) {
        if ( p.verbose ) {
            fprintf(stdout, "Writing attenuation weights to: %s\n", p.saveAttenuationPath);
        }
        FILE *file;
        file = fopen(p.saveAttenuationPath, "wb");
        
        for ( int i=0; i<p.fftParams.paddedT; i=i+1 ) {
            fprintf(file,"%.10f\t%.10f\n", p.fftParams.frequencyBins[i], p.fftParams.binAttenuation[i]);
        }
        
        fclose(file);
    }
    
    p.parametersValid = TRUE;
    return p;
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
    
    struct rsFFTFilterParams fftParams = rsFFTFilterInit(T, T, sampling_rate, 0.01, 0.04, RSFFTFILTER_CUTOFF, 0.0, FALSE, FALSE);
    rsFFTFilter(fftParams, data);
    
    for (int i=0; i<T; i=i+1) {
        fprintf(stderr, "%.10f\n", data[i]);
    }
}