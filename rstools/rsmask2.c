//
//  rsfillholes.c
//  rstools
//
//  Created by André Hoffmann on 8/5/13.
//
//

#include <stdio.h>
#include "rsmathutils.h"

void rsMaskPrintHelp() {
    printf(
        "basic usage:  rsmask2 -input <volume> -output <volume> -mask <volume>\n"
        "\n"
    );
    
    printf(
        "options:\n"
    );
    
    printf(
        "   -input <volume>        : a 3D/4D volume to which the mask will be applied\n"
    );
    
    printf(
        "   -output <volume>       : the volume in which the result will be saved\n"
    );
    
    printf(
        "   -mask <mask>           : a mask specifying the ROI\n"
    );
    
    printf(
        "   -threads <int>         : number of threads used for processing\n"
    );
    
    printf(
        "   -v[erbose]             : show debug information\n"
        "\n"
    );
}

int main(int argc, char * argv[]) {
    
    FSLIO *fslio;
    FSLIO *fslioOutput;
	void *buffer;
    void *outputBuffer;
	size_t buffsize;
	
	char *inputpath = NULL;
	char *outputpath = NULL;
	char *maskpath = NULL;
	
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
        
    BOOL verbose = FALSE;
    int threads = 1;
	
	int ac;
    
	if( argc < 2 ) {
        rsMaskPrintHelp();
        return 1;
    }
    
	/* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if( ! strncmp(argv[ac], "-h", 2) ) {
			rsMaskPrintHelp();
            return 1;
		} else if ( ! strcmp(argv[ac], "-input") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -input\n");
				return 1;
			}
			inputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-m", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m[ask]\n");
				return 1;
			}
			maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-v", 2) ) {
			verbose = TRUE;
		} else if ( ! strcmp(argv[ac], "-output") ) {
            if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -output\n");
				return 1;
			}
			outputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-threads") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -threads\n");
           		return 1;
           	}
           	threads = atoi(argv[ac]);  /* no string copy, just pointer assignment */
        } else {
			fprintf(stderr, "\nError, unrecognized command %s\n", argv[ac]);
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
	
	if ( maskpath == NULL ) {
		fprintf(stderr, "A binary mask must be specified!\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Input file:  %s\n", inputpath);
        fprintf(stdout, "Mask file:   %s\n", maskpath);
        fprintf(stdout, "Output file: %s\n", outputpath);
    }
    
    // Load file
    
    /* open input file */
    fslio = FslOpen(inputpath, "rb");
    if (fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", inputpath);
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
    
    /* prepare centrality file */    
    fslioOutput = FslOpen(outputpath, "wb");
    
    if (fslioOutput == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", outputpath);
        return 1;
    }
    
    FslCloneHeader(fslioOutput, fslio);
    FslSetDim(fslioOutput, xDim, yDim, zDim, vDim);
    FslSetDimensionality(fslioOutput, 4);
    FslSetDataType(fslioOutput, pixtype);
    FslWriteHeader(fslioOutput);
    
    /* load mask */
    unsigned long nPoints = 0L;
    double ***mask = d3matrix(zDim, yDim, xDim);
    Point3D *maskPoints = rsReadMask(maskpath, xDim, yDim, zDim, &nPoints, NULL, fslio, mask);
    if ( maskPoints == NULL) {
        fprintf(stderr, "\nError: Mask invalid.\n");
        FslClose(fslio);
        FslClose(fslioOutput);
        return 1;
    }
    free(mask);
    
    // Prepare buffer
    buffsize       = (size_t)xDim*(size_t)yDim*(size_t)zDim*(size_t)vDim*(size_t)dt/(size_t)8;
    buffer         = malloc(buffsize);
    outputBuffer   = malloc(buffsize);
        
    if (buffer == NULL || outputBuffer == NULL) {
        fprintf(stdout, "Not enough free memory :-(\n");
        return 1;
    }
    
    FslReadVolumes(fslio, buffer, vDim);
    
    rsResetBufferToValue(dt, outputBuffer, slope, inter, xDim, yDim, zDim, vDim, threads, log(-1)); 
    
    long p;
    
    #pragma omp parallel num_threads(threads) private(p) shared(buffer,fslio)
    {
        #pragma omp for schedule(guided, 1)
        for (p = 0L; p<nPoints; p=p+1L) {

            const Point3D point = maskPoints[p];
            double *pointValues = (double*)malloc(sizeof(double)*vDim);
            
            rsExtractTimecourseFromBuffer(fslio, pointValues, buffer, slope, inter, point, xDim, yDim, zDim, vDim);
            rsWriteTimecourseToBuffer(fslioOutput, pointValues, outputBuffer, slope, inter, point, xDim, yDim, zDim, vDim);
            
            free(pointValues);
        }
    }
    
    /* write output */
    free(buffer);
    FslWriteVolumes(fslioOutput, outputBuffer, vDim);
    FslClose(fslioOutput);
    FslClose(fslio);
    free(outputBuffer);
    free(fslioOutput);
    free(maskPoints);
}
