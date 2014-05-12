//
//  rsfillholes.c
//  rstools
//
//  Created by André Hoffmann on 8/5/13.
//
//

#include <stdio.h>
#include "maths/rsmathutils.h"
#include "rsregression_common.h"

void rsMaskPrintHelp() {
    printf(
	   RSTOOLS_VERSION_LABEL "\n"
        "basic usage:  rsmotionscrubbing -input <volume> -rp <rp.txt> -output <volume> -mask <volume>\n"
        "\n"
    );
    
    printf(
        "options:\n"
    );
    
    printf(
        "   -input <volume>        : the 4D volume to be scrubbed\n"
    );
    
    printf(
        "   -output <volume>       : the volume in which the result will be saved\n"
    );

	printf(
		"   -rp <rp.txt>           : the file containing the realignment parameters\n"
	);
	
	printf(
		"   -dvars <dvars.txt>     : (optional) file where the DVARS values will be\n"
		"                            saved to\n"
	);
    
	printf(
		"   -fd <fd.txt>           : (optional) file to which the framewise displacement\n"
		"                            will be saved to\n"
	);
	
	printf(
		"   -flagged <flagged.txt> : (optional) file to which the indices of all\n"
		"                            flagged frames will be saved to\n"
	);
	
	printf(
		"   -dvarsthreshold <float>: (optional) DVARs threshold\n"
	);
	
	printf(
		"   -fdthreshold <float>   : (optional) framewise displacement threshold\n"
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

extern inline double deg2mm(const double dist, const double deg) {
	return M_PI * dist * (deg/180.0);
}

extern inline int rsMax(const int a, const int b) {
  return a > b ? a : b;
}

extern inline int rsMin(const int a, const int b) {
  return a < b ? a : b;
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
	char *realignmentpath = NULL;
	char *fdpath = NULL;
	char *dvarspath = NULL;
	char *flaggedpath = NULL;
	
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim;
	size_t pixtype;
	short dt;
    float inter = 0.0, slope = 1.0;

	float fdthreshold=1.0, dvarsthreshold=0.05;
        
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
		} else if ( ! strcmp(argv[ac], "-rp") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -rp\n");
				return 1;
			}
			realignmentpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-m", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m[ask]\n");
				return 1;
			}
			maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-dvars") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -dvars\n");
				return 1;
			}
			dvarspath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-fd") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -fd\n");
				return 1;
			}
			fdpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-flagged") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -flagged\n");
				return 1;
			}
			flaggedpath = argv[ac];  /* no string copy, just pointer assignment */
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
        } else if ( ! strcmp(argv[ac], "-dvarsthreshold") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -dvarsthreshold\n");
           		return 1;
           	}
           	dvarsthreshold = atof(argv[ac]);  /* no string copy, just pointer assignment */
        } else if ( ! strcmp(argv[ac], "-fdthreshold") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -fdthreshold\n");
           		return 1;
           	}
           	fdthreshold = atof(argv[ac]);  /* no string copy, just pointer assignment */
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
	
	if ( realignmentpath == NULL ) {
		fprintf(stderr, "A realignment parameter file must be specified!\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Input file:  %s\n", inputpath);
        fprintf(stdout, "Mask file:   %s\n", maskpath);
        fprintf(stdout, "Output file: %s\n", outputpath);
		fprintf(stdout, "DVARS threshold: %.4f\n", dvarsthreshold);
		fprintf(stdout, "Framewise displacement threshold: %.4f\n", fdthreshold);
        fprintf(stdout, "Realignment parameters file: %s\n", realignmentpath);
		
		if ( dvarspath != NULL ) 
			fprintf(stdout, "DVARS file: %s\n", dvarspath);
		
		if ( fdpath != NULL ) 
			fprintf(stdout, "Framewise displacement file: %s\n", fdpath);
		
		if ( flaggedpath != NULL ) 
			fprintf(stdout, "Flagged frames file: %s\n", flaggedpath);
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

	/* load realignment parameters */
	long rpColumns;
	long rpEntries;
    double **rp = rsLoadRegressors(realignmentpath, &rpColumns, &rpEntries, 1.0);
	double fd[rpEntries];
	
	/* compute framewise displacement */
	fd[0] = 0;
	for ( int i=1; i<rpEntries; i=i+1 ) {
		//fprintf(stdout, "%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", rp[1][i], rp[2][i], rp[3][i], rp[4][i], rp[5][i], rp[6][i]);
		const double dx     = fabs(rp[1][i]-rp[1][i-1]);
		const double dy     = fabs(rp[2][i]-rp[2][i-1]);
		const double dz     = fabs(rp[3][i]-rp[3][i-1]);
		const double dpitch = fabs(deg2mm(50, rp[4][i])-deg2mm(50, rp[4][i-1]));
		const double droll  = fabs(deg2mm(50, rp[5][i])-deg2mm(50, rp[5][i-1]));
		const double dyaw   = fabs(deg2mm(50, rp[6][i])-deg2mm(50, rp[6][i-1]));
		fd[i] = dx + dy + dz + dpitch + droll + dyaw;
		//fprintf(stdout, "%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", dx, dy, dz, dpitch, droll, dyaw);
		//fprintf(stdout, "%.4f\n", fd[i]);
	}
	
	/* determine datatype */
	pixtype = FslGetDataType(fslio, &dt);
    
    // Prepare buffer
    buffsize       = rsGetBufferSize(xDim, yDim, zDim, vDim, dt);
    buffer         = malloc(buffsize);
        
    if (buffer == NULL) {
        fprintf(stderr, "Not enough free memory :-(\n");
        return 1;
    }

    FslReadVolumes(fslio, buffer, vDim);
    
	/* load mask */
    unsigned long nPoints = 0L;
    double ***mask = d3matrix(zDim, yDim, xDim);
    Point3D *maskPoints = ReadMask(maskpath, xDim, yDim, zDim, &nPoints, NULL, fslio, mask);
    if ( maskPoints == NULL) {
        fprintf(stderr, "\nError: Mask invalid.\n");
        FslClose(fslio);
        FslClose(fslioOutput);
        return 1;
    }
	
	/* get value range */
	double min=9999999999999999, max=-9999999999999999;
	for (t=5; t<vDim; t=t+1) {
		double ***data = d3matrix(zDim-1,  yDim-1, xDim-1);
		rsExtractVolumeFromBuffer(fslio, data[0][0], buffer, slope, inter, t, xDim, yDim, zDim);
		
		for (long p = 0L; p<nPoints; p=p+1L) {
			const Point3D point  = maskPoints[p];
			const double intensity = data[point.z][point.y][point.x];
			
			if ( intensity != intensity ) {
				continue;
			}
			
			if ( min > intensity ) {
				min = intensity;
			}
			
			if ( max < intensity ) {
				max = intensity;
			}
		}
		
		free(data[0][0]); free(data[0]); free(data);
	}
	
	fprintf(stdout, "min: %.2f, max: %.2f\n", min, max);

	/* compute dvars */
	double dvars[vDim];
	dvars[0]=0;
    #pragma omp parallel num_threads(threads) private(t) shared(buffer,fslio,dvars,nPoints,xDim,yDim,zDim,max,min,maskPoints)
    {
        #pragma omp for schedule(guided, 1)
		for (t=1; t<vDim; t=t+1) {

			/* Extract a single volume for timepoint t from the buffer */
			double ***dataNow    = d3matrix(zDim-1,  yDim-1, xDim-1);
			double ***dataBefore = d3matrix(zDim-1,  yDim-1, xDim-1);
	        rsExtractVolumeFromBuffer(fslio, dataNow[0][0],    buffer, slope, inter, t,   xDim, yDim, zDim);
	        rsExtractVolumeFromBuffer(fslio, dataBefore[0][0], buffer, slope, inter, t-1, xDim, yDim, zDim);
			
			dvars[t]=0;
			
	        for (long p = 0L; p<nPoints; p=p+1L) {
				const Point3D point  = maskPoints[p];
				const double INow    = (dataNow[point.z][point.y][point.x]    - min)/fabs(max-min);
				const double IBefore = (dataBefore[point.z][point.y][point.x] - min)/fabs(max-min);
				
				//fprintf(stdout, "%d: X%03dY%03dZ%03d: %.3f - %.3f\n", t, point.x, point.y, point.z, IBefore, INow);
				dvars[t]=dvars[t] + pow(INow - IBefore, 2.0);
	        }
			
			//fprintf(stdout, "%d: %.4f\n", t, dvars[t]);
			dvars[t]=sqrt(dvars[t] / (double)nPoints);
			//fprintf(stdout, "%d: %.4f\n", t, dvars[t]);
	
			/* Free up memory */
			free(dataNow[0][0]);    free(dataNow[0]);	 free(dataNow);
			free(dataBefore[0][0]); free(dataBefore[0]); free(dataBefore);
		}
    }

	if ( dvarspath != NULL ) {
		
		FILE *dvarsfile; 
		dvarsfile = fopen(dvarspath,"w");
		
		for ( t=0; t<vDim; t=t+1 ) {
			fprintf(dvarsfile, "%.10f\n", dvars[t]);
		}
		
		fclose(dvarsfile);
	}

	if ( fdpath != NULL ) {
		
		FILE *fdfile; 
		fdfile = fopen(fdpath,"w");
		
		for ( t=0; t<vDim; t=t+1 ) {
			fprintf(fdfile, "%.10f\n", fd[t]);
		}
		
		fclose(fdfile);
	}
    
	/* mark all frames as unflagged */
	BOOL flaggedFrames[vDim];
	
	for ( t=0; t<vDim; t=t+1 ) {
		flaggedFrames[t] = FALSE;
	}
	
	/* determine which frames need to be flagged */
	for ( t=1; t<vDim; t=t+1 ) {
		if ( fd[t] > fdthreshold || dvars[t] > dvarsthreshold ) {
			flaggedFrames[rsMax(t-1,0)]      = TRUE;
			flaggedFrames[t]                 = TRUE;
			flaggedFrames[rsMin(t+1,vDim-1)] = TRUE;
		}
	}
	
	if ( flaggedpath != NULL ) {
		
		FILE *flaggedfile; 
		flaggedfile = fopen(flaggedpath,"w");
		
		for ( t=0; t<vDim; t=t+1 ) {
			if ( flaggedFrames[t] == TRUE ) {
				fprintf(flaggedfile, "%d\n", t);
			}
		}
		
		fclose(flaggedfile);
	}
	
	/* count remaining frames */
	int remainingFrames = 0;
	int framemap[vDim];

	for ( t=0; t<vDim; t=t+1 ) {
		if ( flaggedFrames[t] == FALSE ) {
			framemap[remainingFrames] = t;
			remainingFrames = remainingFrames + 1;
		}
	}

	/* prepare output file */
	 
    fslioOutput = FslOpen(outputpath, "wb");
    
    if (fslioOutput == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", outputpath);
        return 1;
    }
    
    FslCloneHeader(fslioOutput, fslio);
    FslSetDim(fslioOutput, xDim, yDim, zDim, remainingFrames);
    FslSetDimensionality(fslioOutput, 4);
    FslSetDataType(fslioOutput, dt);
	char *callString = rsMergeStringArray(argc, argv);
    rsWriteNiftiHeader(fslioOutput, callString);
	free(callString);

   	/* create resulting file */
	buffsize     = rsGetBufferSize(xDim, yDim, zDim, remainingFrames, dt);
    outputBuffer = malloc(buffsize);
        
    if (outputBuffer == NULL) {
        fprintf(stderr, "Not enough free memory :-(\n");
        return 1;
    }

    #pragma omp parallel num_threads(threads) private(t) shared(flaggedFrames, fslio, buffer, outputBuffer)
    {
        #pragma omp for schedule(guided, 1)
		for ( t=0; t<remainingFrames; t=t+1 ) {
					
			/* Extract a single volume for timepoint t from the buffer */
			double ***data = d3matrix(zDim-1,  yDim-1, xDim-1);
            rsExtractVolumeFromBuffer(fslio, data[0][0], buffer, slope, inter, framemap[t], xDim, yDim, zDim);

			/* Write back to the buffer */
            rsWriteVolumeToBuffer(fslioOutput, data[0][0], outputBuffer, slope, inter, t, xDim, yDim, zDim);

			/* Free up memory */
			free(data[0][0]); free(data[0]); free(data);
		}
	}

    free(buffer);
    FslWriteVolumes(fslioOutput, outputBuffer, remainingFrames);
    FslClose(fslioOutput);
    free(fslioOutput);
    free(outputBuffer);
    free(maskPoints);
}
