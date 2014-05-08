//
//  rsfastcentrality2.c
//  rstools
//
//  Created by André Hoffmann on 8/8/13.
//
//

#include <stdio.h>
#include "rsmathutils.h"

void rsTestPowerIteration();
void rsTestInvErr();

void rsCentralityPrintHelp() {
    printf(
 	    RSTOOLS_VERSION_LABEL "\n\n"
        "basic usage:  rsfastcentrality2 -input <volume> -output <volume> -mask <volume>\n"
        "\n"
    );
    
    printf(
        "options:\n"
    );
    
    printf(
        "   -input <volume>        : a 4D volume for which the centrality will be computed\n"
    );
    
    printf(
        "   -output <volume>       : the volume in which the result will be saved\n"
    );
    
    printf(
        "   -mask <mask>           : a mask specifying the ROI\n"
    );
 
    printf(
        "   -savemask <mask>       : optional path where the rescaled mask specified with -mask\n"
        "                            will be saved. The saved file with have the same dimensions\n"
        "                            as the input volume.\n"
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
    FSLIO *fslioCentrality;
	void *buffer;
	size_t buffsize;
	
	char *inputpath = NULL;
	char *outputpath = NULL;
	char *maskpath = NULL;
    char *savemaskpath = NULL;
	
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
    int threads = 1;
    
    BOOL verbose = FALSE;
	
	int ac;
    
	if( argc < 2 ) {
        rsCentralityPrintHelp();
        return 1;
    }
    
	/* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if( ! strncmp(argv[ac], "-h", 2) ) {
			rsCentralityPrintHelp();
            return 1;
		} else if ( ! strcmp(argv[ac], "-input") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -input\n");
				return 1;
			}
			inputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-m", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m\n");
				return 1;
			}
			maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-savemask") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -savemask\n");
				return 1;
			}
			savemaskpath = argv[ac];  /* no string copy, just pointer assignment */
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
        } else if ( ! strcmp(argv[ac], "-testPowerIteration") ) {
            rsTestPowerIteration();
            return 0;
        } else if ( ! strcmp(argv[ac], "-testInvErr") ) {
            rsTestInvErr();
            return 0;
        } else {
			fprintf(stderr, "\nError, unrecognized command %s\n", argv[ac]);
		}
	}
	
	rsSetThreadsNum(threads);
	
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
   	void *centralityBuffer;
    
    fslioCentrality = FslOpen(outputpath, "wb");
    
    if (fslioCentrality == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", outputpath);
        return 1;
    }
    
    FslCloneHeader(fslioCentrality, fslio);
    FslSetDim(fslioCentrality, xDim, yDim, zDim, 1);
    FslSetDimensionality(fslioCentrality, 4);
    FslSetDataType(fslioCentrality, pixtype);
	char *callString = rsMergeStringArray(argc, argv);
    rsWriteNiftiHeader(fslioCentrality, callString);
	free(callString);
    
    /* load mask */
    unsigned long nPoints = 0L;
    double ***mask = d3matrix(zDim, yDim, xDim);
    Point3D *maskPoints = ReadMask(maskpath, xDim, yDim, zDim, &nPoints, savemaskpath, fslio, mask);
    if ( maskPoints == NULL) {
        fprintf(stderr, "\nError: Mask invalid.\n");
        FslClose(fslio);
        FslClose(fslioCentrality);
        return 1;
    }
        
    // Prepare buffer
    buffsize = (size_t)xDim*(size_t)yDim*(size_t)zDim*(size_t)vDim*(size_t)dt/(size_t)8;
    buffer   = malloc(buffsize);
    
    if (buffer == NULL) {
        fprintf(stdout, "Not enough free memory :-(\n");
        return 1;
    }
    
    if ( verbose ) fprintf(stdout, "Reading input volume..\n");
    FslReadVolumes(fslio, buffer, vDim);
    FslClose(fslio);
    free(fslio);
    
    long n;
	long nNanPoints = 0;
	long *nanPoints = malloc(sizeof(long)*nPoints);
	
	/* search for voxels that are nan */
	if ( verbose ) fprintf(stdout, "Removing voxels that are NaN..\n");
  	
	for (n=0; n<nPoints; n=n+1) {
        const Point3D point = maskPoints[n];
		
        // load signal
        double *signalData = malloc(sizeof(double) * vDim);
        rsExtractTimecourseFromBuffer(fslioCentrality, signalData, buffer, slope, inter, point, xDim, yDim, zDim, vDim);
        
        for ( int t=0; t<vDim; t=t+1 ) {
			if ( signalData[t] != signalData[t] ) {
				nanPoints[nNanPoints] = n;
				nNanPoints = nNanPoints + 1;
		        free(signalData);
				break;
			}
        }
        
        free(signalData);
	}
	
	/* remove voxels that are nan */
	
	long newIndex = 0;
	Point3D *newMaskPoints = malloc(sizeof(Point3D)*(nPoints-nNanPoints));
	
	for ( n=0; n<nPoints; n=n+1 ) {
		
		if ( ! rsVectorContains(nanPoints, nNanPoints, n) ) {
			newMaskPoints[newIndex] = maskPoints[n];
			newIndex = newIndex + 1;
		}
	}
	free(maskPoints);
	maskPoints = newMaskPoints;
	nPoints = nPoints-nNanPoints;
	
	if ( verbose ) fprintf(stdout, "%lu voxels removed..\n", nNanPoints);

    
    // Implementation according to
    // Wink, Alle Meije, et al.
    // "Fast eigenvector centrality mapping of voxel-wise connectivity in functional magnetic resonance imaging: implementation, validation, and interpretation."
    // Brain Connectivity 2.5 (2012): 265-274.
    
    /* create scaled voxel-timeseries matrix */
    if ( verbose ) fprintf(stdout, "De-meaning %lu voxels followed by rescaling them to have unit variance..\n", nPoints);
    
    gsl_matrix *M  = gsl_matrix_alloc(vDim, nPoints);
    
    const double epsilon  = 2.0e-16;
    const double epsilon2 = 2.0e-8;
    const double corrnorm = sqrt(vDim-1);
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(n)
    {
        #pragma omp for schedule(guided)
        for (n=0; n<nPoints; n=n+1) {
            const Point3D point = maskPoints[n];
            
            // load signal
            double *signalData = malloc(sizeof(double) * vDim);
            rsExtractTimecourseFromBuffer(fslioCentrality, signalData, buffer, slope, inter, point, xDim, yDim, zDim, vDim);
            
            // remove mean and convert to unit variance
            const double mean = gsl_stats_mean(&signalData[0], 1, vDim);
            const double std  = gsl_stats_sd_m(&signalData[0], 1, vDim, mean) + epsilon; // epsilon prevents division by 0
            
            for ( int t=0; t<vDim; t=t+1 ) {
                const double value = (signalData[t] - mean) / (std*corrnorm);
                gsl_matrix_set(M, t, n, value);
            }
            
            free(signalData);
        }
    }
    free(buffer);

    // run power iteration method to get the first eigenvector
    if ( verbose ) fprintf(stdout, "Computing first eigenvector..\n");
    const unsigned int maxIterations = 500;
    unsigned int iteration = 0;
    gsl_vector *vOld = gsl_vector_alloc(nPoints);
    gsl_vector *vNew = gsl_vector_alloc(nPoints);
    gsl_vector_set_all(vOld, 0.5); // inital eigenvector
    gsl_vector *tmpT = gsl_vector_alloc(vDim);
    gsl_vector *tmpN = gsl_vector_alloc(nPoints);
    
    gsl_vector *rightPart = gsl_vector_alloc(nPoints);
    
    while ( iteration < maxIterations ) {
        // start with the result of the last iteration
        gsl_vector_memcpy(vNew, vOld);
        
        // compute the left part of the iteration step formula
        const double leftPart = cblas_dzasum(vDim, vNew->data, 1) / 2.0; // sum(v) / 2
        
        // compute the left part of the iteration step formula
        gsl_blas_dgemv(CblasNoTrans, 1.0, M, vNew, 0.0, tmpT); // tmpT = M^T * v
        gsl_blas_dgemv(CblasTrans, 1.0, M, tmpT, 0.0, rightPart); // rightPart = M * (M^T * v)
        gsl_vector_scale(rightPart, 0.5); // // rightPart = M * (M^T * v) / 2
        
        // add the two parts
        gsl_vector_memcpy(vNew, rightPart);
        gsl_vector_add_constant(vNew, leftPart);
        
        // normalize
        const double invNorm = 1.0 / gsl_blas_dnrm2(vNew);
        gsl_vector_scale(vNew, invNorm);
        
        // check if it has converged
        gsl_vector_sub(vOld, vNew);
        const double diffNorm = gsl_blas_dnrm2(vOld);
        gsl_vector_memcpy(vOld, vNew);
        gsl_vector_scale(vOld, epsilon2);
        const double norm = gsl_blas_dnrm2(vOld);
        gsl_vector_memcpy(vOld, vNew);
        
        iteration = iteration + 1;
        if ( verbose ) fprintf(stdout, "Iteration %03d/%03d Norm Difference: %.17f Target: %.17f\n", iteration, maxIterations, diffNorm, norm);
        
        if (norm > diffNorm) {
            break;
        }
    }
    
    gsl_vector_free(rightPart);
    gsl_vector_free(vOld);
    gsl_matrix_free(M);
    
    // Rescale node centrality
    if ( verbose ) fprintf(stdout, "Rescale values to obey a Gaussian normal distribution\n");
    double *vScaled = malloc(sizeof(double)*nPoints);
    rsSpearmannRank(vScaled, vNew->data, nPoints); // rank data
        
    const double mu    = 0.0;               // mean
    const double sigma = 1.0;               // standard deviation
    const double scale = 2.0 / (nPoints+1); // scaling factor transforming values from uniform [1,N] to uniform ]0,2[
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(n)
    {
        #pragma omp for schedule(guided)
        for (n=0; n<nPoints; n=n+1) {
            vScaled[n] = scale*vScaled[n] - 1; // uniform [1,N] to ]0,2[ and then to ]-1,1[
        }
    }
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(n)
    {
        #pragma omp for schedule(guided)
        for (n=0; n<nPoints; n=n+1) {
            vScaled[n] = rsErfInv(vScaled[n]) * sigma * sqrt(2.0) + mu; // uniform ]-1,1[ to N(mu,sigma)
        }
    }
    
    // Prepare buffer
    buffsize = (size_t)xDim*(size_t)yDim*(size_t)zDim*(size_t)dt/(size_t)8;
    buffer   = malloc(buffsize);
    
    // Write back to file
    if ( verbose ) fprintf(stdout, "Writing output..\n");
    
    rsResetBufferToValue(fslioCentrality->niftiptr->datatype, buffer, slope, inter, xDim, yDim, zDim, 1, sqrt(-1.0));
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(n)
    {
        #pragma omp for schedule(guided)
        for (n = 0L; n<nPoints; n=n+1L) {
            const Point3D point      = maskPoints[n];
            const double  centrality = vScaled[n];
            rsWriteTimecourseToBuffer(fslioCentrality, &centrality, buffer, slope, inter, point, xDim, yDim, zDim, 1);
        }
    }
    
    gsl_vector_free(vNew);
    FslWriteVolumes(fslioCentrality, buffer, 1);
    FslClose(fslioCentrality);
    free(fslioCentrality);
    free(buffer);
    free(maskPoints);
    free(vScaled);
}

void rsTestPowerIteration() {
    
    /*
     >> a = rand(5);
     >> b = triu(a) + triu(a,1)';
     b =
     
     0.2760    0.4984    0.7513    0.9593    0.8407
     0.4984    0.9597    0.2551    0.5472    0.2543
     0.7513    0.2551    0.5060    0.1386    0.8143
     0.9593    0.5472    0.1386    0.1493    0.2435
     0.8407    0.2543    0.8143    0.2435    0.9293
    */
    
    double **m = d2matrix(5,5);
    
    // 1st row
    m[0][0] = 0.2760;
    m[0][1] = 0.4984; m[1][0] = 0.4984;
    m[0][2] = 0.7513; m[2][0] = 0.7513;
    m[0][3] = 0.9593; m[3][0] = 0.9593;
    m[0][4] = 0.8407; m[4][0] = 0.8407;
    
    // 2nd row
    m[1][1] = 0.9597;
    m[1][2] = 0.2551; m[2][1] = 0.2551;
    m[1][3] = 0.5472; m[3][1] = 0.5472;
    m[1][4] = 0.2543; m[4][1] = 0.2543;
    
    // 3rd row
    m[2][2] = 0.5060;
    m[2][3] = 0.1386; m[3][2] = 0.1386;
    m[2][4] = 0.8143; m[4][2] = 0.8143;
    
    // 4th row
    m[3][3] = 0.1493;
    m[3][4] = 0.2435; m[4][3] = 0.2435;
    
    // 5th row
    m[4][4] = 0.9293;
    
    fprintf(stdout, "Input Matrix:\n");
    rs_matrix_fprintf(stdout, (const double**)m, 4, 4, "%.4f");
    fprintf(stdout, "\n");
    
    long double *b = rsFirstEigenvector((const double **)m, 4, 100000, 0.000001, TRUE);
    
    fprintf(stdout, "\n");
    fprintf(stdout, "Eigenvector:\n");
    rs_vector_fprintfl(stdout, b, 4, "%.4Lf");
    
    free(m[0]);
    free(m);
    free(b);
}

void rsTestInvErr()
{
    for (double i=-1; i<=1; i=i+0.25) {
        fprintf(stdout, "%.2f -> %.5f\n", i, rsErfInv(i));
    }
}