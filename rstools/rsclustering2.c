//
//  rsclustering.c
//  rstools
//
//  Created by André Hoffmann on 10/1/13.
//
//

#include <stdio.h>
#include "rsmathutils.h"

void rsClusteringPrintHelp() {
    printf(
        "basic usage:  rsclustering2 -input <volume> -output <volume> -mask <volume>\n"
        "\n"
    );
    
    printf(
        "options:\n"
    );
    
    printf(
        "   -input <volume>        : a 4D volume that will be clustered\n"
    );
    
    printf(
        "   -output <volume>       : the volume in which the result will be saved\n"
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
    FSLIO *fslioClustered;
	void *buffer;
	void *clusteringBuffer;
	size_t buffsize;
	
	char *inputpath = NULL;
	char *outputpath = NULL;
	char *maskpath = NULL;
	
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
    
    short kernelsize = 15;
    
    BOOL verbose = FALSE;
    int threads = 1;
	
	int ac;
    
	if( argc < 2 ) {
        rsClusteringPrintHelp();
        return 1;
    }
    
	/* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if( ! strncmp(argv[ac], "-h", 2) ) {
			rsClusteringPrintHelp();
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
	
    if ( verbose ) {
        fprintf(stdout, "Input file:  %s\n", inputpath);
        fprintf(stdout, "Output file: %s\n", outputpath);
        fprintf(stdout, "Mask file:   %s\n", maskpath);
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
    fslioClustered = FslOpen(outputpath, "wb");
    
    if (fslioClustered == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", outputpath);
        return 1;
    }
    
    FslCloneHeader(fslioClustered, fslio);
    FslSetDim(fslioClustered, xDim, yDim, zDim, vDim);
    FslSetDimensionality(fslioClustered, 4);
    FslSetDataType(fslioClustered, pixtype);
    FslWriteHeader(fslioClustered);
    
	// Prepare buffer
    buffsize          = (size_t)xDim*(size_t)yDim*(size_t)zDim*(size_t)vDim*(size_t)dt/(size_t)8;
    buffer   		  = malloc(buffsize);
        
    if (buffer == NULL) {
        fprintf(stdout, "Not enough free memory :-(\n");
        return 1;
    }
    
    FslReadVolumes(fslio, buffer, vDim);

	
	/* load mask */
    unsigned long nPoints = 0L;
    double ***mask = d3matrix(zDim, yDim, xDim);
	
	const double threshold = 0.01;
	void *thresholdPointer = (void*) &threshold;
	
    Point3D *maskPoints = rsReadMaskCustomThreshold(maskpath, xDim, yDim, zDim, &nPoints, NULL, fslio, mask, rsNegativeThreshold, thresholdPointer);
    if ( maskPoints == NULL) {
        fprintf(stderr, "\nError: Mask invalid.\n");
        FslClose(fslio);
		FslClose(fslioClustered);
        return 1;
    }
	free(mask);
    
	long p;
	
	#pragma omp parallel num_threads(threads) private(p) shared(buffer,fslio)
    {
        #pragma omp for schedule(guided, 1)
        for (p = 0L; p<nPoints; p=p+1L) {
            const Point3D point = maskPoints[p];
    		
    		double *signal = malloc(sizeof(double)*vDim);
    		
    		for (t = 0; t<vDim; t=t+1) {
    			signal[t] = log(-1.0);
    		}
   
            /* write out z-scored data to buffer */
            rsWriteTimecourseToBuffer(fslio, &signal[0], buffer, slope, inter, point, xDim, yDim, zDim, vDim);

			free(signal);
    	}
    }

	/* Convert the input to z-scores */
	#pragma omp parallel num_threads(threads) private(y,x,t) shared(buffer,fslio, slope, inter, xDim, yDim, zDim)
    {
        #pragma omp for schedule(guided)
        for (z=0; z<zDim; z=z+1) {
            for (y=0; y<yDim; y=y+1) {
                for (x=0; x<xDim; x=x+1) {
                    
					Point3D point = MakePoint3D(x, y, z);
                    
                    /* read out timecourse from buffer */
                    double *signal = malloc(vDim*sizeof(double));
                    rsExtractTimecourseFromBuffer(fslio, signal, buffer, slope, inter, point, xDim, yDim, zDim, vDim);
                    
					double mean = gsl_stats_mean(signal, 1, vDim);
					double std  = gsl_stats_sd_m(signal, 1, vDim, mean);
					
					for (t = 0; t<vDim; t=t+1) {
						signal[t] = (signal[t] - mean) / std;
					}

                    /* write out z-scored data to buffer */
                    rsWriteTimecourseToBuffer(fslio, signal, buffer, slope, inter, point, xDim, yDim, zDim, vDim);
                    free(signal);
                }
            }
        }
    }

    #pragma omp parallel num_threads(threads) private(t,x,y,z) shared(buffer,fslio,clusteringBuffer,fslioClustered, slope, inter, xDim, yDim, zDim)
    {
        #pragma omp for schedule(guided, 1)
        for (t = 0; t<vDim; t=t+1) {
            
			/* Extract a single volume for timepoint t from the buffer */
			double ***data = d3matrix(zDim-1,  yDim-1, xDim-1);
            rsExtractVolumeFromBuffer(fslio, data[0][0], buffer, slope, inter, t, xDim, yDim, zDim);

			

			/* Write back to the buffer */
            rsWriteVolumeToBuffer(fslioClustered, data[0][0], buffer, slope, inter, t, xDim, yDim, zDim);

			/* Free up memory */
			free(data[0][0]);
			free(data[0]);
			free(data);
        }
    }
    
    /* write output */
    FslClose(fslio);
    free(fslio);
    FslWriteVolumes(fslioClustered, buffer, vDim);
    FslClose(fslioClustered);
    free(fslioClustered);
    free(buffer);
}