//
//  rscorrelation_common.c
//  rstools
//
//  Created by André Hoffmann on 7/9/13.
//
//

#include "rscorrelation_common.h"

void rsCorrelationPrintHelp() {
    printf(
	   RSTOOLS_VERSION_LABEL "\n"
       "rscorrelation: This tool will correlate the timecourse of\n"
       "               of every voxel in the supplied volume with\n"
       "               a timecourse that is supplied via standard\n"
       "               input."
       "\n"
    );
    
    printf(
       "basic usage:  rscorrelation -input <volume> -mask <volume> -output <volume> [-verbose]\n"
       "\n"
    );
    
    printf(
       "options:\n"
    );
    
    printf(
       "   -help                  : show this help\n"
    );
    
    printf(
       "   -input <volume>        : the volume for which the correlation of the\n"
       "                            timecourse for each voxel is computed\n"
    );
    
    printf(
       "   -output <volume>       : the volume in which the correlation values will\n"
       "                            be saved in\n"
    );
    
    printf(
       "   -mask <mask>           : a mask specifying the ROI for improved performance\n"
    );
    
	printf(
	   "   -montecarlo <n> <m>    : repeats the computation of the correlation n times and uses\n"
	   "                            m randomly drawn timepoint in each run. eventually the average\n"
	   "                            is being saved. (using it enforces the conversion to z-scores)\n"
	);
    
    printf(
       "   -conversion <mode>     : <mode> specifies what is stored in the output file, it can\n"
       "                            take the follwing values: \n"
       "                             * none: the correlation coefficients will be stored without\n"
       "                                     converting them\n"
       "                             * z:    the correlation coefficients will be converted to\n"
       "                                     z-values before being stored(default)\n"
       "                             * t:    the correlation coefficients will be converted to\n"
       "                                     values of the T-statistic\n"
    );

    printf(
       "   -threads <int>         : number of threads used for processing\n"
    );
    
    printf(
       "   -comment <txt>         : adds a comment about the origin of thereference timecourse\n"
       "                            to the nifti header of the correlation map\n"
    );

    printf(
       "   -delay <int>           : delay the regressor by the defined volumes(or TR)\n"
    );
    
    printf(
       "   -v[erbose]             : show debug information\n"
       "\n"
    );
}

struct rsCorrelationParameters rsCorrelationInitParameters() {
    struct rsCorrelationParameters p;
    
    p.inputpath             = NULL;
    p.maskpath              = NULL;
    p.outputpath            = NULL;
    p.savemaskpath          = NULL;
	p.commentpath           = NULL;
	p.comment               = NULL;
    p.xDim                  = 0;
    p.yDim                  = 0;
    p.zDim                  = 0;
    p.vDim                  = 0;
    p.delay                 = 0;
    p.pixtype               = 0;
    p.dt                    = 4;
    p.inter                 = 0.0;
    p.slope                 = 1.0;
    p.verbose               = FALSE;
    p.fslio                 = NULL;
    p.fslioCorrelation      = NULL;
    p.correlation           = NULL;
    p.parametersValid       = FALSE;
    p.mask                  = NULL;
    p.nRegressorValues      = 0;
    p.regressor             = NULL;
    p.threads               = 1;
    p.conversionMode        = RSTOOLS_CORRELATION_CONVERSION_Z;
	p.monteCarloRepetitions = 0;
	p.monteCarloSampleSize  = 0;


    return p;
}

struct rsCorrelationParameters rsCorrelationLoadParams(int argc, char * argv[]) {

    struct rsCorrelationParameters p = rsCorrelationInitParameters();
    
    int ac;
        
	/* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if ( ! strcmp(argv[ac], "-input") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -input\n");
				return p;
			}
			p.inputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-output") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -output\n");
				return p;
			}
			p.outputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-montecarlo") ) {
  			if( (++ac)+1 >= argc ) {
           		fprintf(stderr, "** missing arguments for -montecarlo\n");
           		return p;
           	}
           	p.monteCarloRepetitions = atoi(argv[ac]);
			++ac;
           	p.monteCarloSampleSize  = atoi(argv[ac]);
	    } else if ( ! strncmp(argv[ac], "-m", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m\n");
				return p;
			}
			p.maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-s", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -savemask\n");
				return p;
			}
			p.savemaskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-comment") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -comment\n");
				return p;
			}
			p.commentpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-conversion") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -conversion\n");
				return p;
			}
            
            if ( ! strcmp(argv[ac], "none") ) {
                p.conversionMode = RSTOOLS_CORRELATION_CONVERSION_NONE;
            } else if ( ! strcmp(argv[ac], "z") ) {
                p.conversionMode = RSTOOLS_CORRELATION_CONVERSION_Z;
            } else if ( ! strcmp(argv[ac], "t") ) {
                p.conversionMode = RSTOOLS_CORRELATION_CONVERSION_T;
            } else {
                fprintf(stderr, "** invalid value for -conversion\n");
				return p;
            }
		} else if ( ! strcmp(argv[ac], "-threads") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -threads\n");
           		return p;
           	}
           	p.threads = atoi(argv[ac]);  /* no string copy, just pointer assignment */
        } else if ( ! strcmp(argv[ac], "-delay") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -delay\n");
           		return p;
           	}
           	p.delay = atoi(argv[ac]);
        } else if ( ! strncmp(argv[ac], "-v", 2) ) {
			p.verbose = TRUE;
		} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
			return p;
		}
	}
	
	if ( p.inputpath == NULL ) {
		fprintf(stderr, "No input volume specified!\n");
		return p;
	}
	
	if ( p.outputpath == NULL ) {
		fprintf(stderr, "No output volume specified!\n");
		return p;
	}
	
	if ( p.commentpath != NULL ) {
		p.comment = rsReadCommentFile(p.commentpath);
	}
	
	char *callString = rsMergeStringArray(argc, argv);
	char *comment;
	if ( p.comment == NULL ) {
		comment = callString;
	} else {
		char *separator = "\nReference Timecourse Info:\n";
		comment = malloc(sizeof(char)*(strlen(callString)+strlen(separator)+strlen(p.comment)+1));
		sprintf(&comment[0], "%s%s%s", callString, separator, p.comment);
		free(callString);
	}
	
    /* Read seed timecourse from stdin */
    p.regressor = rsReadRegressorFromStandardInput(&p.nRegressorValues);
    
    if ( p.verbose ) {
        fprintf(stdout, "Input file: %s\n", p.inputpath);
        fprintf(stdout, "Mask file: %s\n", p.maskpath);
        fprintf(stdout, "Seed length: %u\n", p.nRegressorValues);
        fprintf(stdout, "Seed delay: %u\n", p.delay);
    }
    
    /* open input file */
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
	p.dt = FslGetDataType(p.fslio, &p.pixtype);
    
    /* prepare correlation file */
   	void *correlationBuffer;
    
    p.fslioCorrelation = FslOpen(p.outputpath, "wb");
    
    if (p.fslioCorrelation == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", p.outputpath);
        return p;
    }
    
    FslCloneHeader(p.fslioCorrelation, p.fslio);
    FslSetDim(p.fslioCorrelation, p.xDim, p.yDim, p.zDim, 1);
    FslSetDimensionality(p.fslioCorrelation, 4);
    FslSetDataType(p.fslioCorrelation, p.pixtype);
    rsWriteNiftiHeader(p.fslioCorrelation, comment);
    
    /* prepare buffer */
    p.correlation = d3matrix(p.zDim, p.yDim, p.xDim);
    
    /* load mask */
    p.mask = NULL;
    if ( p.maskpath != NULL ) {
        unsigned long nPoints = 0L;
        p.mask = d3matrix(p.zDim, p.yDim, p.xDim);
        Point3D *maskPoints = ReadMask(p.maskpath, p.xDim, p.yDim, p.zDim, &nPoints, p.savemaskpath, p.fslio, p.mask);
        if ( maskPoints == NULL) {
            fprintf(stderr, "\nError: Mask invalid.\n");
            FslClose(p.fslio);
            FslClose(p.fslioCorrelation);
            return p;
        }
        free(maskPoints);
    }
    
    p.parametersValid = TRUE;
    
    return p;
}

void rsCorrelationWriteCorrelationFile(struct rsCorrelationParameters* p) {
    
    /* Write correlation file */
    if ((*p).verbose) fprintf(stdout, "Writing correlation file\n");
    size_t buffsize = (size_t)((size_t)(*p).xDim*(size_t)(*p).yDim*(size_t)(*p).zDim*(size_t)(*p).dt/(size_t)8);
    void *correlationBuffer;
    correlationBuffer = malloc(buffsize);
    
    convertScaledDoubleToBuffer(
        (*p).fslioCorrelation->niftiptr->datatype,
        correlationBuffer,
        &(*p).correlation[0][0][0],
        (*p).fslioCorrelation->niftiptr->scl_slope,
        (*p).fslioCorrelation->niftiptr->scl_inter,
        (*p).xDim,
        (*p).yDim,
        (*p).zDim,
        TRUE
    );
    
    FslWriteVolumes((*p).fslioCorrelation, correlationBuffer, 1);
    
    free(correlationBuffer);
}

void rsCorrelationFree(struct rsCorrelationParameters* p) {
    FslClose((*p).fslioCorrelation);
    FslClose((*p).fslio);
    free((*p).fslioCorrelation);
    free((*p).fslio);
}

double *rsReadRegressorFromStandardInput(unsigned int *nValues) {
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    double value;
    unsigned int nBuffer = 10;
    *nValues = 0;
    double *regressor = (double*)malloc(nBuffer * sizeof(double));
    
    while ((read = getline(&line, &len, stdin)) != -1) {
        value = atof(line);
        regressor[*nValues] = value;
        *nValues = *nValues + 1;
        
        // Check if we're running out of memory and extend the array if necessary
        if ( *nValues + 1 >= nBuffer ) {
            nBuffer = nBuffer + 10;
            double* tmpRegressor = realloc(regressor, nBuffer * sizeof(double));
            if (tmpRegressor) {
                regressor = tmpRegressor;
            } else {
                fprintf(stderr, "Could not allocate enough memory to save the regressor from stdin.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    
    if (line) free(line);
    
    return regressor;
}