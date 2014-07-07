//
//  rscentrality2.c
//  rstools
//
//  Created by André Hoffmann on 7/23/13.
//
//

#include <stdio.h>
#include "src/maths/rsmathutils.h"

void rsTestPowerIteration();
void rsCentralityPrintHelp() {
    printf(
  	    RSTOOLS_VERSION_LABEL "\n"
        "basic usage:  rscentrality -input <volume> -output <volume> -mask <volume>\n"
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
        "   -correlation <mode>    : specifies how the correlation coefficients are converted before\n"
        "                            the eigenvalue computation. <mode> can be on of the following:\n"
        "                            pos, neg, abs(default), scaled\n"
    );
    
    printf(
        "   -savemask <mask>       : optional path where the rescaled mask specified with -mask\n"
        "                            will be saved. The saved file with have the same dimensions\n"
        "                            as the input volume.\n"
    );
    
    printf(
        "   -savesimilarity <file> : optional path to a file where the similarity matrix will be \n"
        "                            saved to.\n"
    );
    
    printf(
        "   -similarity <file>     : optional path to use an existing similarity matrix instead of\n"
        "                            creating it file where the similarity matrix will be \n"
    );
    
    
    printf(
        "   -precise               : uses long double instead of double for internal computations\n"
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
    char *savesimilaritypath = NULL;
    char *similaritypath = NULL;
	
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim;
	size_t pixtype;
	short dt;
    float inter = 0.0, slope = 1.0;
    int threads = 1;
    int correlationMode = RSMATRIXCONVERSION_ABSOLUTE;
    
    BOOL verbose = FALSE;
    BOOL precise = FALSE;
	
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
		} else if ( ! strcmp(argv[ac], "-precise") ) {
			precise = TRUE;
		} else if ( ! strcmp(argv[ac], "-output") ) {
            if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -output\n");
				return 1;
			}
			outputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-savesimilarity") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -savesimilarity\n");
				return 1;
			}
			savesimilaritypath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-similarity") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -similarity\n");
				return 1;
			}
			similaritypath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-threads") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -threads\n");
           		return 1;
           	}
           	threads = atoi(argv[ac]);  /* no string copy, just pointer assignment */
        } else if ( ! strcmp(argv[ac], "-correlation") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -correlation\n");
           		return 1;
           	} else if ( ! strcmp(argv[ac], "abs") ) {
                correlationMode = RSMATRIXCONVERSION_ABSOLUTE;
            } else if ( ! strcmp(argv[ac], "pos") ) {
                correlationMode = RSMATRIXCONVERSION_POSITIVE;
            } else if ( ! strcmp(argv[ac], "neg") ) {
                correlationMode = RSMATRIXCONVERSION_NEGATIVE;
            } else if ( ! strcmp(argv[ac], "scaled") ) {
                correlationMode = RSMATRIXCONVERSION_SCALED;
            } else {
                fprintf(stderr, "** invalid argument for -correlation(must be abs, pos or neg)\n");
           		return 1;
            }
        } else if ( ! strcmp(argv[ac], "-test") ) {
            rsTestPowerIteration();
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
        if (savesimilaritypath)
            fprintf(stdout, "Similarity matrix output file: %s\n", savesimilaritypath);
        if (similaritypath)
            fprintf(stdout, "Similarity matrix file: %s\n", similaritypath);
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
   	void *centralityBuffer;
    
    fslioCentrality = FslOpen(outputpath, "wb");
    
    if (fslioCentrality == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", outputpath);
        return 1;
    }
    
    FslCloneHeader(fslioCentrality, fslio);
    FslSetDim(fslioCentrality, xDim, yDim, zDim, 1);
    FslSetDimensionality(fslioCentrality, 4);
    FslSetDataType(fslioCentrality, dt);
	char *callString = rsMergeStringArray(argc, argv);
    rsWriteNiftiHeader(fslioCentrality, callString);
	free(callString);
    
    /* load mask */
    unsigned long nPoints = 0L;
    double ***mask = d3matrix(zDim-1, yDim-1, xDim-1);
    Point3D *maskPoints = rsReadMask(maskpath, xDim, yDim, zDim, &nPoints, savemaskpath, fslio, mask);
    if ( maskPoints == NULL) {
        fprintf(stderr, "\nError: Mask invalid.\n");
        FslClose(fslio);
        FslClose(fslioCentrality);
        return 1;
    }
    
    double **similarity = d2matrix(nPoints-1, nPoints-1);
    
    if (similaritypath) {
        /* If similarity matrix exists load it */
        
        if ( rsLoadMatrix(similaritypath, similarity, nPoints, nPoints) != TRUE ) {
            fprintf(stderr, "Error while saving similarity matrix '%s'\n", savesimilaritypath);
            FslClose(fslio);
            FslClose(fslioCentrality);
            return 1;
        }
        
        // Prepare buffer
        buffsize = rsGetBufferSize(xDim, yDim, zDim, 1, dt);
        buffer   = rsMalloc(buffsize);
        
        if (buffer == NULL) {
            fprintf(stdout, "Not enough free memory :-(\n");
            return 1;
        }
        
    } else {
        /* If similarity matrix doesn't exist yet, create it */
        
        // Prepare buffer
        buffsize = rsGetBufferSize(xDim, yDim, zDim, vDim, dt);
        buffer   = malloc(buffsize);
        
        if (buffer == NULL) {
            fprintf(stdout, "Not enough free memory :-(\n");
            return 1;
        }
        
        FslReadVolumes(fslio, buffer, vDim);
        
        /* create similarity matrix */
        if ( verbose ) fprintf(stdout, "Creating the similarity matrix for %lu voxels..\n", nPoints);
        double *signal1;
        double *signal2;
        long p1, p2;
        long processedPoints = 0L;
        
        #pragma omp parallel num_threads(threads) private(p1,p2,signal1,signal2) shared(buffer,fslio,similarity)
        {
            #pragma omp for schedule(guided, 1)
            for (p1 = 0L; p1<nPoints; p1++) {
                for (p2 = p1; p2<nPoints; p2++) {
                    
                    const Point3D *point1 = &maskPoints[p1];
                    const Point3D *point2 = &maskPoints[p2];
                    
                    /* read out timecourses from buffer */
                    signal1 = malloc(vDim*sizeof(double));
                    signal2 = malloc(vDim*sizeof(double));
                    rsExtractTimecourseFromBuffer(dt, signal1, buffer, slope, inter, point1, xDim, yDim, zDim, vDim);
                    rsExtractTimecourseFromBuffer(dt, signal2, buffer, slope, inter, point2, xDim, yDim, zDim, vDim);
                    
                    /* compute correlation */
                    double correlation;
                    
                    if ( precise ) {
                        correlation = rsZCorrelation(signal1, signal2, vDim);
                    } else {
                        correlation = rsFastZCorrelation(signal1, signal2, vDim);
                    }
                    similarity[p1][p2] = correlation;
                    similarity[p2][p1] = correlation;
                    
                    /* cleanup */
                    free(signal1);
                    free(signal2);
                }
                
                /* show progress */
                #pragma omp atomic
                processedPoints += 1L;
                
                if (verbose && processedPoints > 0 && processedPoints % (long)(nPoints / 10) == 0) {
                    fprintf(stdout, "..%.0f%%\n", ceil((float)processedPoints*100.0 / (float)nPoints));
                }
            }
        }
    }
    
    FslClose(fslio);
    free(fslio);
    
    /* save similarity matrix */
    if ( savesimilaritypath != NULL ) {
        if ( verbose ) fprintf(stdout, "Saving the similarity matrix..\n");
        
        if ( rsSaveMatrix(savesimilaritypath, (const double **)similarity, nPoints, nPoints) != TRUE ) {
            fprintf(stderr, "Error while saving similarity matrix '%s'\n", savesimilaritypath);
        }
    }
    
    /* convert correlation coefficients in the similarity matrix */
    if ( verbose ) {
        fprintf(
            stdout,
            "Converting correlation coefficients to %s values..\n",
            (correlationMode == RSMATRIXCONVERSION_ABSOLUTE ? "absolute" : (correlationMode == RSMATRIXCONVERSION_POSITIVE ? "strictly positive" : (correlationMode == RSMATRIXCONVERSION_NEGATIVE ? "strictly negative" : "rescaled" )))
        );
    }
    rsMatrixConversion(similarity, nPoints, nPoints, correlationMode);
    
    /* compute first eigenvector */
    if ( verbose ) fprintf(stdout, "Computing eigenvector centrality..\n");
    
    long double *eigenvector = rsFirstEigenvector((const double**)similarity, nPoints, 10000, 0.000001, verbose);
    free(similarity[0]);
    free(similarity);
    
    /* write back to file */
    if ( verbose ) fprintf(stdout, "Writing output..\n");
    
    rsResetBufferToValue(fslioCentrality->niftiptr->datatype, buffer, slope, inter, xDim, yDim, zDim, 1, sqrt(-1.0));
    
    for (unsigned long p = 0L; p<nPoints; p=p+1L) {
        Point3D *point     = &maskPoints[p];
        double centrality = fabs(eigenvector[p]);
        
        rsWriteTimecourseToBuffer(dt, &centrality, buffer, slope, inter, point, xDim, yDim, zDim, 1);
    }
    
    FslWriteVolumes(fslioCentrality, buffer, 1);
    FslClose(fslioCentrality);
    free(fslioCentrality);
    free(buffer);
    free(eigenvector);
    free(maskPoints);
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
    
    double **m = d2matrix(4,4);
    
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
