//
//  rsregression_common.c
//  rstools
//
//  Created by André Hoffmann on 6/27/13.
//
//

#include "rsregression_common.h"

void rsRegressionPrintHelp() {
   printf(
      "basic usage:  rsregression -input <volume> -regressors <txtFile> [-residuals <volume> | -fitted <volume> | -betas <volume> | -mask <volume>]\n"
      "\n"
   );
    
   printf(
      "options:\n"
   );
 
   printf(
      "   -input <volume>      : the volume to be regressed\n"
   );
    
   printf(
      "   -residuals <volume>  : the volume in which the residuals will be saved\n"
   );
    
   printf(
      "   -fitted <volume>     : the volume in which the fitted volumes will be saved\n"
   );
   
   printf(
      "   -betas <volume>      : the volume to be regressed\n"
   );
    
   printf(
      "   -mask <mask>         : a mask specifying the ROI for improved performance\n"
   );
    
   printf(
      "   -savemask <mask>       : optional path where the rescaled mask specified with -mask\n"
      "                            will be saved. The saved file with have the same dimensions\n"
      "                            as the input volume.\n"
   );
    
   printf(
      "   -regressors <txt>      : a tabbed/spaced textfile containing the regressors with the\n"
      "                            different regressors in the columns and time course in the\n"
      "                            rows. Decimal numbers may be formatted like this: 1.23e+45\n"
   );
    
   printf(
      "   -f1 <double>           : (optional) the lower frequency of the bandpass\n"
      "                            filter\n"
   );

   printf(
      "   -f2 <double>           : (optional) the upper frequency of the bandpass filter\n"
      "                            filter\n"
   );

   printf(
      "   -samplingrate <double> : (optional) only required if bandpass filtering is\n"
      "                            requested.\n"
   );
    
   printf(
      "   -threads <int>         : (rsregression2 only) number of threads used for processing\n"
   );
   
   printf(
      "   -v[erbose]             : show debug information\n"
      "\n"
   );
}

struct rsRegressionParameters rsRegressionInitParameters() {
    struct rsRegressionParameters p;
    
    p.inputpath            = NULL;
    p.maskpath             = NULL;
    p.regressorspath       = NULL;
    p.savemaskpath         = NULL;
    p.saveBetasPath        = NULL;
    p.saveResidualsPath    = NULL;
    p.saveFittedPath       = NULL;
    p.xDim                 = 0;
    p.yDim                 = 0;
    p.zDim                 = 0;
    p.vDim                 = 0;
    p.pixtype              = 0;
    p.dt                   = 4;
    p.inter                = 0.0;
    p.slope                = 1.0;
    p.f1                   = -1.0;
    p.f2                   = -1.0;
    p.sampling_rate        = -1.0;
    p.verbose              = FALSE;
    p.filterActive         = FALSE;
    p.nRegressors          = 0;
    p.nRegressorValues     = 0;
    p.regressors           = NULL;
    p.fslio                = NULL;
    p.parametersValid      = FALSE;
    p.mask                 = NULL;
    p.nyquist_frequency    = 0.0;
    p.bin_width            = 0.0;
    p.nFrequencyBinsLow    = 0;
    p.nFrequencyBinsHigh   = 0;
    p.nFrequencyBins       = 0;
    p.nFrequencyRegressors = 0;
    p.frequencyBins        = 0;
    p.nAllRegressors       = 0;
    p.allRegressors        = 0;
    p.threads              = 1;
    
    return p;
}

struct rsRegressionParameters rsRegressionLoadParams(int argc, char * argv[]) {

    struct rsRegressionParameters p = rsRegressionInitParameters();
    
    /* parse parameters */
	for( int ac = 1; ac < argc; ac++ ) {
		if ( ! strcmp(argv[ac], "-input") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -input\n");
				return p;
			}
			p.inputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-betas") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -betas\n");
				return p;
			}
			p.saveBetasPath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-residuals") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -residuals\n");
				return p;
			}
			p.saveResidualsPath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-fitted") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -fitted\n");
				return p;
			}
			p.saveFittedPath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-regressors") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -regressors\n");
				return p;
			}
			p.regressorspath = argv[ac];  /* no string copy, just pointer assignment */
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
		} else if ( ! strcmp(argv[ac], "-f1") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -f1\n");
				return p;
			}
			p.f1 = atof(argv[ac]);  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-threads") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -threads\n");
				return p;
			}
			p.threads = atoi(argv[ac]);  /* no string copy, just pointer assignment */
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
		} else if ( ! strncmp(argv[ac], "-v", 2) ) {
			p.verbose = TRUE;
		} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
            return p;
		}
	}
	
	if ( p.inputpath == NULL ) {
		fprintf(stderr, "No input volume specified(-input)!\n");
		return p;
	}
	
	if ( p.regressorspath == NULL ) {
		fprintf(stderr, "A file containing the regressors must be specified(-regressors)!\n");
		return p;
	}
    
    p.filterActive = p.sampling_rate >= 0.0 || p.f1 >= 0.0 || p.f2 >= 0.0;
    
    if ( p.verbose ) {
        fprintf(stdout, "Input file: %s\n", p.inputpath);
        fprintf(stdout, "Regressors file: %s\n", p.regressorspath);
        fprintf(stdout, "Mask file: %s\n", p.maskpath);
        fprintf(stdout, "Residuals file: %s\n", p.saveResidualsPath);
        fprintf(stdout, "Fitted file: %s\n", p.saveFittedPath);
        fprintf(stdout, "Betas file: %s\n", p.saveBetasPath);
        
        if ( p.filterActive ) {
            fprintf(stdout, "Bandpass filter active!\n");
            fprintf(stdout, "Sampling rate: %.4f\n", p.sampling_rate);
            fprintf(stdout, "Bandpass: (%.4f-%.4fHz)\n", p.f1, p.f2);
        }
    }
    
    if ( p.filterActive && (p.sampling_rate < 0.0 || p.f1 < 0.0 || p.f2 < 0.0) ) {
        fprintf(stderr, "Bandpass filter parameters are not complete. (-samplingrate, -f1, -f2)!\n");
		return p;
    }
    
    // load regressors
    p.regressors = rsLoadRegressors(p.regressorspath, &p.nRegressors, &p.nRegressorValues, 1.0);
    
    if ( p.verbose ) {
        fprintf(stdout, "Regressors: %ld, Samples: %ld\n", p.nRegressors, p.nRegressorValues);
    }
    
    p.fslio = FslOpen(p.inputpath, "rb");
    if (p.fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", p.inputpath);
        return p;
    }
    
	/* determine dimensions */
	FslGetDim(p.fslio, &p.xDim, &p.yDim, &p.zDim, &p.vDim);
    
    if ( p.verbose ) {
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", p.xDim, p.yDim, p.zDim, p.vDim);
    }
    
    if (p.fslio->niftiptr->scl_slope != 0) {
        p.slope = p.fslio->niftiptr->scl_slope;
        p.inter = p.fslio->niftiptr->scl_inter;
    }
	
	/* determine datatype and initalize buffer */
	p.dt       = FslGetDataType(p.fslio, &p.pixtype);
    p.wordsize = p.fslio->niftiptr->nbyper;
    
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
    
    if ( p.filterActive ) {
    
        // Prepare bandpass filter
        p.nyquist_frequency = 1.0 / (2.0 * p.sampling_rate);
        p.bin_width         = p.nyquist_frequency / (p.vDim/2);
        
        // Compute the number of frequency regressors that will be added to the existing regressors
        p.nFrequencyBinsLow    = (int)floor(p.f1 / p.bin_width);                         // number of frequency bins before the lower cutoff frequency
        p.nFrequencyBinsHigh   = (int)floor((p.nyquist_frequency - p.f2) / p.bin_width); // number of frequency bins after the highter cutoff frequency
        p.nFrequencyBins       = p.nFrequencyBinsLow + p.nFrequencyBinsHigh;             // number of frequency bins in total
        p.nFrequencyRegressors = p.nFrequencyBins * 2;                                   // number of frequency regressors (both sine and cosine)
        
        // Compute the frequencies for the bins
        p.frequencyBins = malloc(p.nFrequencyBins * sizeof(double));
        fprintf(stdout, "Bin width: %.4f\n", p.bin_width);
        for ( int i=0; i<p.nFrequencyBins; i=i+1 ) {
            if ( (i < p.nFrequencyBinsLow) ) {
                p.frequencyBins[i] = (i + 1) * p.bin_width;
            } else {
                p.frequencyBins[i] = (i - p.nFrequencyBinsLow + 1) * p.bin_width + p.f2;
            }
            
            fprintf(stdout, "Frequency bin %d(%s): %.4f\n", i, (i < p.nFrequencyBinsLow ? "low" : "high"), p.frequencyBins[i]);
        }
    }
    
    // Add frequency regressors to the other regressors if requested
    p.nAllRegressors = p.filterActive ? p.nRegressors + 1 + p.nFrequencyRegressors : p.nRegressors + 1;
    
    if ( p.filterActive ) {
        p.allRegressors = d2matrix(p.nAllRegressors, p.vDim);
        
        for (int i=0; i<p.nAllRegressors; i=i+1) {
            
            if ( i < p.nRegressors ) {
                for (int t=0; t<p.vDim; t=t+1) {
                    p.allRegressors[i][t] = p.regressors[i][t];
                }
            } else if ( i < (p.nRegressors + p.nFrequencyBins) ) {
                const int j = i - p.nRegressors;
                for (int t=0; t<p.vDim; t=t+1) {
                    p.allRegressors[i][t] = rsSampleSineWave(p.sampling_rate, p.frequencyBins[j], t);
                }
            } else {
                const int j = i - p.nRegressors - p.nFrequencyBins;
                for (int t=0; t<p.vDim; t=t+1) {
                    p.allRegressors[i][t] = rsSampleCosineWave(p.sampling_rate, p.frequencyBins[j], t);
                }
            }
        }
    } else {
        p.allRegressors = p.regressors;
    }
    
    p.parametersValid = TRUE;
    return p;
}

/*
 * Reads in a single line from a file and returns it.
 */
BOOL rsReadline(FILE *f, char *line, int *length) {
    *length = 0;
    int c;
    
    while(TRUE) {
        c = fgetc(f);
        
        if (c == '\n' || c == '\r' || c == EOF) {
            break;
        }
        
        line[*length] = (char)c;
        *length = *length+1;
    }
    
    *(*line+length) = '\0';
    
    return c!=EOF;
}

/*
 * Reads in a line from the regressor file and returns the
 * tab- or space-separated values as an ar array of doubles.
 */
double *rsParseRegressorLine(char *line, long *nRegressors) {
    double *regressors;
    char delimiter[] = " \t";
    char *ptr;
    *nRegressors = 0L;
    BOOL endsWithSeparator = TRUE;
    
    size_t lineLength = strlen(line);
    
    // if it doesn't end with a tab or space we need to add one for strtok to work
    if (line[lineLength-1] != '\t' || line[lineLength-1] != ' ') {
        lineLength = lineLength + 1;
        endsWithSeparator = FALSE;
    }
    
    char lineCpy[lineLength];
    strcpy(lineCpy, line);
    
    if ( !endsWithSeparator ) {
        lineCpy[lineLength-1] = '\t';
    }
    
    ptr = strtok(lineCpy, delimiter);
    while (ptr != NULL) {
        ptr = strtok(NULL, delimiter);
        *nRegressors = *nRegressors+1;
    }
    
    regressors = malloc(sizeof(double)*(*nRegressors));
    
    strcpy(lineCpy, line);
    
    if ( !endsWithSeparator ) {
        lineCpy[lineLength-1] = '\t';
    }
    
    long n=0L;
    ptr = strtok(lineCpy, delimiter);
    while(TRUE) {
        if (ptr == NULL) {
            break;
        }
        
        regressors[n] = atof(ptr);
        n = n+1L;
        ptr = strtok(NULL, delimiter);
    }
    
    return regressors;
}

/*
 * Loads a regressor file where the regressors are in the columns
 * and the samples in the rows. Regressor values should be separated
 * by spaces or tabs.
 * Returns a 2D matrix in the following form: double[regressor][time]
 */
double **rsLoadRegressors(char *path, long *nRegressors, long *nValues, double constantFactor) {
    FILE *f = fopen(path, "r");
    
    if (f == NULL) {
        fprintf(stderr, "Error: Regressors could not be read.\n");
        return NULL;
    }
    
    /* Read regressor file to get the number of regressors and samples */
    char *line = malloc(sizeof(char)*1000);
    int length = 0;
    *nRegressors = -1L;
    *nValues = 0L;
    
    rewind(f);
    
    long n = 0L;
    while( rsReadline(f, line, &length) ) {
        if ( length < 1 ) {
            continue;
        }
        
        double *regressors = rsParseRegressorLine(line, &n);
        if (*nRegressors < 0) {
            *nRegressors = n;
        }
        free(regressors);
        *nValues = *nValues+1L;
    }
    
    /* Initialize result matrix */
    rewind(f);
    double **result = d2matrix(*nRegressors+1, *nValues);
    
    /* Fill with constant regressor */
    for ( long v=0L; v<n; v = v+1L ) {
        result[0][v] = constantFactor;
    }
    
    /* Save regressors in result matrix */
    long v=0L;
    while( rsReadline(f, line, &length) ) {
        if ( length < 1 ) {
            continue;
        }
        
        double *regressors = rsParseRegressorLine(line, &n);
        
        for( long l=0L; l<n; l=l+1L ) {
            result[l+1L][v] = regressors[l];
        }
        free(regressors);
        v=v+1L;
    };
    
    free(line);
    fclose(f);
    
    return result;
}