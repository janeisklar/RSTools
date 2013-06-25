//
//  rscorrelation.c
//  rstools
//
//  Created by André Hoffmann on 6/7/13.
//
//
#include <stdio.h>
#include <nifti1.h>
#include <fslio.h>
#include "rsniftiutils.h"
#include "rsmathutils.h"

double *readRegressorFromStandardInput(unsigned int *nValues);

int show_help( void )
{
    printf(
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
       "   -v[erbose]             : show debug information\n"
       "\n"
    );
    
    return 0;
}

int main(int argc, char * argv[]) {

    FSLIO *fslio;
	void *buffer;
	size_t buffsize;
	
	char *inputpath    = NULL;
	char *outputpath   = NULL;
	char *maskpath     = NULL;
    char *savemaskpath = NULL;
	
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
    
    BOOL verbose = FALSE;
	
	int ac;
    
	if( argc < 2 ) return show_help();   /* typing '-help' is sooo much work */
    
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
		} else if ( ! strcmp(argv[ac], "-output") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -output\n");
				return 1;
			}
			outputpath = argv[ac];  /* no string copy, just pointer assignment */
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
		fprintf(stderr, "No input volume specified!\n");
		return 1;
	}
	
	if ( outputpath == NULL ) {
		fprintf(stderr, "No output volume specified!\n");
		return 1;
	}
	
    /* Read seed timecourse from stdin */
    unsigned int nValues;
    double *regressor = readRegressorFromStandardInput(&nValues);
    
    if ( verbose ) {
        fprintf(stdout, "Input file: %s\n", inputpath);
        fprintf(stdout, "Mask file: %s\n", maskpath);
        fprintf(stdout, "Seed length: %u\n", nValues);
    }
    
    /* open input file */
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
    
    /* prepare correlation file */
    FSLIO *fslioCorrelation;
   	void *correlationBuffer;
    
    fslioCorrelation = FslOpen(outputpath, "wb");
    
    if (fslioCorrelation == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n",outputpath);
        return 1;
    }
    
    FslCloneHeader(fslioCorrelation, fslio);
    FslSetDim(fslioCorrelation, xDim, yDim, zDim, 1);
    FslSetDimensionality(fslioCorrelation, 4);
    FslSetDataType(fslioCorrelation, pixtype);
    FslWriteHeader(fslioCorrelation);
    
    /* prepare buffer */
    double ***correlation = d3matrix(zDim,yDim,xDim);
    buffsize = vDim*dt/8;
    buffer = malloc(buffsize);
    
    /* load mask */
    double ***mask = NULL;
    if ( maskpath != NULL ) {
        unsigned long nPoints = 0L;
        mask = d3matrix(zDim, yDim, xDim);
        Point3D *maskPoints = ReadMask(maskpath, xDim, yDim, zDim, &nPoints, savemaskpath, fslio, mask);
        if ( maskPoints == NULL) {
            fprintf(stderr, "\nError: Mask invalid.\n");
            FslClose(fslio);
            FslClose(fslioCorrelation);
            return 1;
        }
        free(maskPoints);
    }

    /* Iterate over all voxels for which the correlation is to be computed */
    double signal[vDim];
    for (short z=0; z<zDim; z=z+1) {
        if (verbose) fprintf(stdout, "Computing slice Z%03hd/%03hd\n", z+1, zDim);
        for (short y=0; y<yDim; y=y+1) {
            for (short x=0; x<xDim; x=x+1) {
                                
                /* If it's not in the mask skip it to improve the performance */
                if (mask != NULL && mask[z][y][x] < 0.1) {
                    
                    /* set the value in the correlation file to NaN so that it is skipped in later processing steps */
                    correlation[z][y][x] = log(-1.0);

                    continue;
                }
                
                /* read out timecourse */
                FslReadTimeSeries(fslio, buffer, x, y, z, vDim);
                convertBufferToScaledDouble(signal, buffer, (long)vDim, slope, inter, fslio->niftiptr->datatype);
                
                /* compute correlation */
                correlation[z][y][x] = rsZCorrelation(signal, regressor, (size_t)nValues);
            }
        }
    }
    
    /* Write correlation file */
    if (verbose) fprintf(stdout, "Writing correlation file\n");
    buffsize = (size_t)((size_t)xDim*(size_t)yDim*(size_t)zDim*(size_t)dt/(size_t)8);
    correlationBuffer = malloc(buffsize);
    
    convertScaledDoubleToBuffer(
        fslioCorrelation->niftiptr->datatype,
        correlationBuffer,
        &correlation[0][0][0],
        fslioCorrelation->niftiptr->scl_slope,
        fslioCorrelation->niftiptr->scl_inter,
        xDim,
        yDim,
        zDim,
        TRUE
    );
    
    FslWriteVolumes(fslioCorrelation, correlationBuffer, 1);
    
    /* Close everything */
    FslClose(fslioCorrelation);
    FslClose(fslio);
    free(fslioCorrelation);
    free(fslio);
    
	return 0;
}

double *readRegressorFromStandardInput(unsigned int *nValues) {
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