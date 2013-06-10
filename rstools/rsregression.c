/*******************************************************************
 *
 * rstimecourse.c
 *
 * Extracts the time course from a 4D-Nifti by supplying either a voxel or a binary mask
 *
 * Usage: rstimecourse [-m <mask> [-a <algorithm>]] [-p <X> <Y> <Z>] -input <volume>
 *
 * 
 * André Hoffmann
 *******************************************************************/


#include <stdio.h>
#include <strings.h>
#include <regex.h>

#include <nifti1.h>
#include <fslio.h>
#include "rsniftiutils.h"
#include "rsmathutils.h"

BOOL rsReadline(FILE *f, char *line, int *length);
double **rsLoadRegressors(char *path, long *nRegressors, long *nValues, double constantFactor);

int show_help( void )
{
   printf(
      "rsregression: Given a 4D-Nifti and a txt file with regressors(columns),\n"
      "              this tool will perform a multiple linear regression on it.\n"
      "\n"
   );
    
   printf(
      "basic usage:  rsregression -input <volume> -regressors <txtFile> [-residuals <volume> | -fitted <volume> | -betas <volume> | -mask <volume>]\n"
      "\n"
   );
    
   printf(
      "options:\n"
   );

   printf(
      "   -help                : show this help\n"
   );
 
   printf(
      "   -input <volume>      : the volume to be regressed\n"
   );
    
   printf(
      "   -residuals <volume>  : the volume in which the residuals will be saved\n"
   );
    
   printf(
      "   -fitted <volume>     : the volume in which the fitted volumes will be saved\n"
      "                          !!! This function has not been implemented yet !!!\n"
   );
   
   printf(
      "   -betas <volume>      : the volume to be regressed\n"
   );
    
   printf(
      "   -mask <mask>         : a mask specifying the ROI for improved performance\n"
   );
    
   printf(
      "   -savemask <mask>     : optional path where the rescaled mask specified with -mask\n"
      "                          will be saved. The saved file with have the same dimensions\n"
      "                          as the input volume.\n"
   );
    
   printf(
      "   -regressors <txt>    : a tabbed/spaced textfile containing the regressors with the\n"
      "                          different regressors in the columns and time course in the\n"
      "                          rows. Decimal numbers may be formatted like this: 1.23e+45\n"
   );
   
   printf(
      "   -v[erbose]           : show debug information\n"
      "\n"
   );
    
   return 0;
}

int main(int argc, char * argv[])
{
    FSLIO *fslio;
	void *buffer;
	size_t buffsize;
	
	char *inputpath = NULL;
	char *maskpath = NULL;
	char *regressorspath = NULL;
    char *savemaskpath = NULL;
    char *saveBetasPath = NULL;
    char *saveResidualsPath = NULL;
    char *saveFittedPath = NULL;
	
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
    
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
		} else if ( ! strcmp(argv[ac], "-betas") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -betas\n");
				return 1;
			}
			saveBetasPath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-residuals") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -residuals\n");
				return 1;
			}
			saveResidualsPath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-fitted") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -fitted\n");
				return 1;
			}
			saveFittedPath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-regressors") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -regressors\n");
				return 1;
			}
			regressorspath = argv[ac];  /* no string copy, just pointer assignment */
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
		} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
		}
	}
	
	if ( inputpath == NULL ) {
		fprintf(stderr, "No input volume specified(-input)!\n");
		return 1;
	}
	
    /*
	if ( maskpath == NULL ) {
		fprintf(stderr, "A binary mask must be specified(-mask)!\n");
		return 1;
	}*/
	
	if ( regressorspath == NULL ) {
		fprintf(stderr, "A file containing the regressors must be specified(-regressors)!\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Input file: %s\n", inputpath);
        fprintf(stdout, "Mask file: %s\n", maskpath);
        fprintf(stdout, "Residuals file: %s\n", saveResidualsPath);
        fprintf(stdout, "Fitted file: %s\n", saveFittedPath);
        fprintf(stdout, "Betas file: %s\n", saveBetasPath);
    }
    
    // load regressors
    long nRegressors = 0L;
    long nRegressorValues = 0L;

    double **regressors;
    regressors = rsLoadRegressors(regressorspath, &nRegressors, &nRegressorValues, 1.0);
    
    if ( verbose ) {
        fprintf(stdout, "Regressors: %ld, Samples: %ld\n", nRegressors, nRegressorValues);
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

    /* prepare residuals file */
    FSLIO *fslioResiduals;
   	void *residualsBuffer;
    
    if ( saveResidualsPath != NULL ) {

        fslioResiduals = FslOpen(saveResidualsPath, "wb");
        
        if (fslioResiduals == NULL) {
            fprintf(stderr, "\nError, could not read header info for %s.\n",saveResidualsPath);
            return 1;
        }
        
        FslCloneHeader(fslioResiduals, fslio);
        FslSetDim(fslioResiduals, xDim, yDim, zDim, vDim);
        FslSetDimensionality(fslioResiduals, 4);
        FslSetDataType(fslioResiduals, pixtype);
        FslWriteHeader(fslioResiduals);
        
        // prepare buffer
        buffsize = (size_t)((size_t)vDim*(size_t)dt/(size_t)8);
        residualsBuffer = malloc(buffsize);

	if (verbose) fprintf(stdout, "residualsbuffer: %lu\n", buffsize);
    }
    
    /* prepare betas file */
    FSLIO *fslioBetas;
   	void *betasBuffer;
    
    if ( saveBetasPath != NULL ) {
        
        fslioBetas = FslOpen(saveBetasPath, "wb");
        
        if (fslioBetas == NULL) {
            fprintf(stderr, "\nError, could not read header info for %s.\n",saveBetasPath);
            return 1;
        }
        
        FslCloneHeader(fslioBetas, fslio);
        FslSetDim(fslioBetas, xDim, yDim, zDim, (nRegressors+1L));
        FslSetDimensionality(fslioBetas, 4);
        FslSetDataType(fslioBetas, pixtype);
        FslWriteHeader(fslioBetas);
        
        // prepare buffer
        buffsize = (size_t)((size_t)xDim*(size_t)yDim*(size_t)zDim*(size_t)(nRegressors+1)*(size_t)dt/(size_t)8);
        betasBuffer = malloc(buffsize);
    }
    
    FSLIO *fslioFitted;
 
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
    buffsize = vDim*dt/8;
    buffer = malloc(buffsize);
    double signal[vDim];
    double betas[nRegressors+1L];
    double residuals[vDim];
    double fitted[vDim];

    /* Prepare empty timecourse */
    int emptyBufferLength = vDim > nRegressors ? vDim : nRegressors+1;
    void *emptybuffer     = malloc((size_t)emptyBufferLength*(size_t)dt/(size_t)8);
    
    double v[vDim];
    for (short t=0; t<emptyBufferLength; t=t+1) {
        v[t] = 0.0;
    }
    convertScaledDoubleToBuffer(fslioResiduals->niftiptr->datatype, emptybuffer, v, slope, inter, emptyBufferLength, 1, 1, FALSE);
    
    /* Iterate over all voxels that are to be regressed */
    BOOL regressionInitalized = FALSE;
    for (short z=0; z<zDim; z=z+1) {
        fprintf(stdout, "Regressing slice Z%03hd/%03hd\n", z+1, zDim);
        for (short y=0; y<yDim; y=y+1) {
            for (short x=0; x<xDim; x=x+1) {
                
                if (verbose) fprintf(stdout, "(%03hd,%03hd,%03hd)\n", x, y, z);
                
                /* If it's not in the mask skip it to improve the performance */
                if (mask != NULL && mask[z][y][x] < 0.1) {
                    
                    /* set the value in the residuals to 0 so that the nifti isn't empty */
                    if ( fslioResiduals != NULL ) {
                        rsWriteTimeSeries(fslioResiduals, emptybuffer, x, y, z, vDim);
                    }

                    /* set the value in the betas to 0 so that the nifti isn't empty */
                    if ( fslioBetas != NULL ) {
                        rsWriteTimeSeries(fslioBetas, emptybuffer, x, y, z, nRegressors+1L);
                    }
                    
                    continue;
                }
                
                /* read out timecourse */
                FslReadTimeSeries(fslio, buffer, x, y, z, vDim);
                convertBufferToScaledDouble(signal, buffer, (long)vDim, slope, inter, fslio->niftiptr->datatype);
                
                /* run the regression */
                rsLinearRegression(
                    (int)vDim,
                    signal,
                    (int)nRegressors+1,
                    regressors,
                    betas,
                    residuals,
                    fitted,
                    verbose
                );
                
                /* write out residuals if desired */
                if ( saveResidualsPath != NULL ) {
                    convertScaledDoubleToBuffer(fslioResiduals->niftiptr->datatype, residualsBuffer, residuals, slope, inter, vDim, 1, 1, FALSE);
                    rsWriteTimeSeries(fslioResiduals, residualsBuffer, x, y, z, vDim);
                }
                
                /* write out betas if desired */
                if ( saveBetasPath != NULL ) {
                    convertScaledDoubleToBuffer(fslioBetas->niftiptr->datatype, betasBuffer, betas, slope, inter, nRegressors+1L, 1, 1, FALSE);
                    rsWriteTimeSeries(fslioBetas, betasBuffer, x, y, z, nRegressors+1L);
                }
            }
        }
    }
    
    if ( saveResidualsPath != NULL ) {
        FslClose(fslioResiduals);
        free(fslioResiduals);
        free(residualsBuffer);    }
    
    if ( saveBetasPath != NULL ) {
        FslClose(fslioBetas);
        free(fslioBetas);
        free(betasBuffer);
    }
    
    if ( maskpath != NULL ) {
        free(mask);
    }
    
    free(buffer);
    FslClose(fslio);
    free(fslio);
    free(regressors);
    free(emptybuffer);
    
	return 0;
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
    
    char lineCpy[strlen(line)];
    strcpy(lineCpy, line);
    
    ptr = strtok(lineCpy, delimiter);
    while (ptr != NULL) {
        ptr = strtok(NULL, delimiter);
        *nRegressors = *nRegressors+1;
    }
    
    regressors = malloc(sizeof(double)*(*nRegressors));
    
    strcpy(lineCpy, line);
    
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

//double rsParse(char *value) {
//
//    regex_t regex;
//    int reti;
//    regmatch_t match[2];
//    double result = 0.0;
//    char msgbuf[100];
//
//    printf("*Current value: %s.\n", value);
//
//    /* Compile regular expression */
//    reti = regcomp(&regex, "^([-+]?[\\d]*\\.?[\\d]*)e([-+]?[\\d]*)", REG_ICASE | REG_EXTENDED);
//
//    if( reti ){
//        fprintf(stderr, "Could not compile regex\n");
//        exit(1);
//    }
//
//    /* Execute regular expression */
//    reti = regexec(&regex, value, 2, match, 0);
//    if( reti == REG_NOMATCH ) {
//        regfree(&regex);
//        return atof(value);
//    } else {
//        regerror(reti, &regex, msgbuf, sizeof(msgbuf));
//        fprintf(stderr, "Regex match failed: %s\n", msgbuf);
//        exit(1);
//    }
//
//    printf("Match #1: %.*s\n", (int)match[0].rm_eo - (int)match[0].rm_so, value + (int)match[0].rm_so);
//    printf("Match #2: %.*s\n", (int)match[1].rm_eo - (int)match[1].rm_so, value + (int)match[1].rm_so);
//
//    /* Free compiled regular expression if you want to use the regex_t again */
//    regfree(&regex);
//}

//void testRegression()
//{
//    int nSamples        = 10;
//    double signal[]     = {1.9, 0.6, 0.5, 0.8, -0.4, -0.9, -0.7, -0.1, -1.7, -0.2};
//    int nRegressors     = 2;
//    /*double regressors[] = {-2.0, -1.0, -0.8, -0.3, 0.0, 0.5, 0.6, 0.7, 1.0, 1.2}; */
//    double **regressors = d2matrix(nRegressors, nSamples);
//    double r1[]         = {-2.0, -1.0, -0.8, -0.3, 0.0, 0.5, 0.6, 0.7, 1.0, 1.2};
//    double r2[]         = {0.3, 0.1, 0.9, 0.1, 1.3, -1.3, -0.7,  0.4, -1.2, 0.8};
//    regressors[0]       = r1;
//    regressors[1]       = r2;
//    double betas[]      = {0.0, 0.0};
//    double residuals[]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//    double fitted[]     = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//
//    rsLinearRegression(nSamples, signal, nRegressors, regressors, betas, residuals, fitted, TRUE);
//
//    for ( int i = 0; i < nSamples; i++) {
//        printf("Residual %d: %.4f\n", i+1, residuals[i]);
//    }
//
//    for ( int i = 0; i < nRegressors+1; i++) {
//        printf("Beta %d: %.4f\n", i+1, betas[i]);
//    }
//
//    for ( int i = 0; i < nSamples; i++) {
//        printf("Fitted %d: %.4f\n", i+1, fitted[i]);
//    }
//}
