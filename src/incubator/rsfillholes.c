//
//  rsfillholes.c
//  rstools
//
//  Created by André Hoffmann on 8/5/13.
//
//

#include <stdio.h>
#include "src/maths/rsmathutils.h"

void rsFillHolesPrintHelp() {
    printf(
 	    RSTOOLS_VERSION_LABEL "\n\n"
        "basic usage:  rsfillholes -input <volume> -output <volume> -mask <volume>\n"
        "\n"
    );
    
    printf(
        "options:\n"
    );
    
    printf(
        "   -input <volume>        : a 3D volume in which the holes will be filled\n"
    );
    
    printf(
        "   -output <volume>       : the volume in which the result will be saved\n"
    );
    
    printf(
        "   -mask <mask>           : a mask specifying the ROI\n"
    );
    
    printf(
        "   -kernel <size>         : specifies how many neighbouring voxels(cube width) should be\n"
        "                            considered when filling in empty values. defaults to 15.\n"
    );
    
    
    printf(
        "   -threads <int>         : number of threads used for processing\n"
    );
    
    printf(
        "   -v[erbose]             : show debug information\n"
        "\n"
    );
}

static inline int max(int a, int b) {
    return a > b ? a : b;
}

static inline int min(int a, int b) {
    return a < b ? a : b;
}

int main(int argc, char * argv[]) {
    
    FSLIO *fslio;
    FSLIO *fslioFilled;
	void *buffer;
	size_t buffsize;
	
	char *inputpath = NULL;
	char *outputpath = NULL;
	char *maskpath = NULL;
	
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim;
	size_t pixtype;
	short dt;
    float inter = 0.0, slope = 1.0;
    
    short kernelsize = 15;
    
    BOOL verbose = FALSE;
    int threads = 1;
	
	int ac;
    
	if( argc < 2 ) {
        rsFillHolesPrintHelp();
        return 1;
    }
    
	/* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if( ! strncmp(argv[ac], "-h", 2) ) {
			rsFillHolesPrintHelp();
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
		} else if ( ! strcmp(argv[ac], "-kernel") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -kernel\n");
				return 1;
			}
			kernelsize = atoi(argv[ac]);
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
        fprintf(stdout, "Kernel size: %d\n", kernelsize);
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
	pixtype = FslGetDataType(fslio, &dt);
    
    /* prepare centrality file */    
    fslioFilled = FslOpen(outputpath, "wb");
    
    if (fslioFilled == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", outputpath);
        return 1;
    }
    
    FslCloneHeader(fslioFilled, fslio);
    FslSetDim(fslioFilled, xDim, yDim, zDim, 1);
    FslSetDimensionality(fslioFilled, 4);
    FslSetDataType(fslioFilled, dt);
	char *callString = rsMergeStringArray(argc, argv);
    rsWriteNiftiHeader(fslioFilled, callString);
	free(callString);
    
    /* load mask */
    unsigned long nPoints = 0L;
    double ***mask = d3matrix(zDim-1, yDim-1, xDim-1);
    Point3D *maskPoints = rsReadMask(maskpath, xDim, yDim, zDim, &nPoints, NULL, fslio, mask);
    if ( maskPoints == NULL) {
        fprintf(stderr, "\nError: Mask invalid.\n");
        FslClose(fslio);
        FslClose(fslioFilled);
        return 1;
    }
    
    // Prepare buffer
    buffsize = rsGetBufferSize(xDim, yDim, zDim, vDim, dt);
    buffer   = malloc(buffsize);
        
    if (buffer == NULL) {
        fprintf(stdout, "Not enough free memory :-(\n");
        return 1;
    }
    
    FslReadVolumes(fslio, buffer, 1);
    
    const short kerneldist = (short)ceil(kernelsize  / 2.0);
    long p;
    
    #pragma omp parallel num_threads(threads) private(p) shared(buffer,fslio)
    {
        #pragma omp for schedule(guided, 1)
        for (p = 0L; p<nPoints; p++) {

            const Point3D *point = &maskPoints[p];
            short kernellength  = 0;
            double pointValue   = 0;
            double norm         = 0;
            double average      = 0;
            
            rsExtractTimecourseFromBuffer(dt, &pointValue, buffer, slope, inter, point, xDim, yDim, zDim, 1);
            
            if ( pointValue == pointValue && pointValue != 0 ) {
                continue;
            }
            
            for (int x=max(point->x-kerneldist,0); x<min(point->x+kerneldist, yDim-1); x=x+1) {
                for (int y=max(point->y-kerneldist,0); y<min(point->y+kerneldist, yDim-1); y=y+1) {
                    for (int z=max(point->z-kerneldist,0); z<min(point->z+kerneldist, zDim-1); z=z+1) {
                        
                        if ( x==point->x && y==point->y && z==point->z ) {
                            continue;
                        }
                        
                        double kernelPointValue = 0;
                        rsExtractTimecourseFromBuffer(dt, &kernelPointValue, buffer, slope, inter, rsMakePoint3D(x, y, z), xDim, yDim, zDim, 1);
                        
                        if ( kernelPointValue != kernelPointValue || kernelPointValue == 0 ) {
                            continue;
                        }
                        
                        const double weight = 3 * kernelsize - abs(point->x-x) - abs(point->y-y) - abs(point->z-z);
                        average      = average + kernelPointValue * weight;
                        norm         = norm + weight;
                        kernellength = kernellength + 1;
                    }
                }
            }
            
            if (kernellength < 1) {
                if(verbose) fprintf(stderr, "Couldn't fill hole for X%03dY%03dZ%03d. Kernel too small?\n", point->x, point->y, point->z);
                continue;
            }
            
            average = average / norm;
            
            rsWriteTimecourseToBuffer(dt, &average, buffer, slope, inter, point, xDim, yDim, zDim, 1);
        }
    }
    
    /* write output */
    
    FslWriteVolumes(fslioFilled, buffer, 1);
    FslClose(fslioFilled);
    free(fslioFilled);
    free(buffer);
    free(maskPoints);
}
