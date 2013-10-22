//
//  rsfillholes.c
//  rstools
//
//  Created by André Hoffmann on 8/5/13.
//
//

#include <stdio.h>
#include "rsmathutils.h"
#include "rsregression_common.h"

void rsMaskPrintHelp() {
    printf(
        "basic usage:  rsmotionscrubbing2 -input <volume> -rp <rp.txt> -output <volume> -mask <volume>\n"
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
	
	if ( realignmentpath == NULL ) {
		fprintf(stderr, "A realignment parameter file must be specified!\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Input file:  %s\n", inputpath);
        fprintf(stdout, "Mask file:   %s\n", maskpath);
        fprintf(stdout, "Output file: %s\n", outputpath);
        fprintf(stdout, "Realignment parameters file: %s\n", realignmentpath);
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
	dt = FslGetDataType(fslio, &pixtype);
    
    /* prepare output file */    
//    fslioOutput = FslOpen(outputpath, "wb");
//    
//    if (fslioOutput == NULL) {
//        fprintf(stderr, "\nError, could not read header info for %s.\n", outputpath);
//        return 1;
//    }
//    
//    FslCloneHeader(fslioOutput, fslio);
//    FslSetDim(fslioOutput, xDim, yDim, zDim, vDim);
//    FslSetDimensionality(fslioOutput, 4);
//    FslSetDataType(fslioOutput, pixtype);
//    FslWriteHeader(fslioOutput);
    
    // Prepare buffer
    buffsize       = (size_t)xDim*(size_t)yDim*(size_t)zDim*(size_t)vDim*(size_t)dt/(size_t)8;
    buffer         = malloc(buffsize);
        
    if (buffer == NULL) {
        fprintf(stderr, "Not enough free memory :-(\n");
        return 1;
    }

    FslReadVolumes(fslio, buffer, vDim);
    
	double signal[vDim];
	double ***dataNow = d3matrix(zDim-1,  yDim-1, xDim-1);
	rsExtractTimecourseFromBuffer(fslio, signal, buffer, slope, inter, MakePoint3D(50, 30, 20), xDim, yDim, zDim, vDim);
    rsExtractVolumeFromBuffer(fslio, dataNow[0][0], buffer, slope, inter, t, xDim, yDim, zDim);
	fprintf(stdout, "%d: %.5f\n", t, dataNow[50][30][20]);
	fprintf(stdout, "%d: %.5f\n", t, signal[t]);
	free(dataNow[0][0]);    free(dataNow[0]);	 free(dataNow);


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

	/* compute dvars */
	double dvars[vDim];
	dvars[0]=0;
    //#pragma omp parallel num_threads(threads) private(t) shared(buffer,fslio,dvars)
    //{
    //    #pragma omp for schedule(guided, 1)
//		for (t=1; t<2/*vDim*/; t=t+1) {
//
//			/* Extract a single volume for timepoint t from the buffer */
//			double ***dataNow    = d3matrix(zDim-1,  yDim-1, xDim-1);
//			double ***dataBefore = d3matrix(zDim-1,  yDim-1, xDim-1);
//	        rsExtractVolumeFromBuffer(fslio, dataNow[0][0],    buffer, slope, inter, t,   xDim, yDim, zDim);
//	        rsExtractVolumeFromBuffer(fslio, dataBefore[0][0], buffer, slope, inter, t-1, xDim, yDim, zDim);
//			
//			fprintf(stdout, "%d: %.5f\n", t, dataNow[50][30][20]);
//			
//			dvars[t]=0;
//			
//	        for (long p = 0L; p<nPoints; p=p+1L) {
//				const Point3D point = maskPoints[p];
//				//fprintf(stdout, "%d: X%03dY%03dZ%03d: %.3f - %.3f\n", t, point.x, point.y, point.z, dataNow[point.z][point.y][point.x], dataBefore[point.z][point.y][point.x]);
//				dvars[t]=dvars[t] + pow(dataNow[point.z][point.y][point.x] - dataBefore[point.z][point.y][point.x], 2.0);
//	        }
//			
//			//fprintf(stdout, "%d: %.4f\n", t, dvars[t]);
//			dvars[t]=sqrt(dvars[t]);
//	
//			/* Free up memory */
//			free(dataNow[0][0]);    free(dataNow[0]);	 free(dataNow);
//			free(dataBefore[0][0]); free(dataBefore[0]); free(dataBefore);
//		}
    //}

	/*for ( t=0; t<vDim; t=t+1 ) {
		fprintf(stdout, "%d: %.4f\n", t, dvars[t]);
	}*/

	
    
    /* write output */
    free(buffer);
    /*FslWriteVolumes(fslioOutput, outputBuffer, vDim);
    FslClose(fslioOutput);
    free(outputBuffer);
    free(fslioOutput);*/
    free(maskPoints);
}
