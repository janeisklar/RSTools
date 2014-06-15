#include "rstimecourse_common.h"

void rsTimecourseInit(rsTimecourseParameters* p)
{
	p->parametersValid = FALSE;

	// check if the required arguments have been provided

	if ( p->inputpath == NULL ) {
		fprintf(stderr, "No input volume specified!\n");
		return;
	}

	if ( p->maskpath == NULL && p->point == NULL ) {
		fprintf(stderr, "Either a binary mask or a voxel coordinate must be specified!\n");
		return;
	}

	if ( p->mask2path == NULL && p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_CSP ) {
		fprintf(stderr, "A second binary mask must be supplied to run CSP! (use --mask2)\n");
		return;
	}

    if ( p->verbose ) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
		if ( p->outputpath != NULL ) {
        	fprintf(stdout, "Output file: %s\n", p->outputpath);
		}
        fprintf(stdout, "Mask file: %s\n", p->maskpath);
		if ( p->mask2path != NULL ) {
        	fprintf(stdout, "Mask2 file: %s\n", p->mask2path);
		}
		if ( p->point != NULL ) {
	    	fprintf(stdout, "Voxel: %d %d %d\n", p->point->x, p->point->y, p->point->z);
		} else if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_MEAN ) {
        	fprintf(stdout, "Algorithm: Mean\n");
		} else if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_STDDEV ) {
			fprintf(stdout, "Algorithm: Standard Deviation\n");
		} else if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_SPCA || p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_TPCA || p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_CSP ) {
			if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_SPCA ) {
				fprintf(stdout, "Algorithm: sPCA\n");
			} else if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_TPCA ) {
				fprintf(stdout, "Algorithm: tPCA\n");
			} else if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_CSP ) {
				fprintf(stdout, "Algorithm: CTP\n");
			}
			if ( p->spatialmappath != NULL ) {
				fprintf(stdout, "Spatial map file: %s\n", p->spatialmappath);
			}
			if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_SPCA || p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_TPCA ) {
				fprintf(stdout, "Variance to be retained: %.2f\n", p->minVariance);
			}
			if ( p->nComponents > 0 ) {
				fprintf(stdout, "Number of components to be retained: %d\n", p->nComponents);
			}
		}
    }

	// if we are not working with a single voxel only read the whole file
	BOOL fileNeedsToBeRead = p->point == NULL;

    // open input file
	p->input = rsOpenNiftiFile(p->inputpath, fileNeedsToBeRead ? RSNIFTI_OPEN_READ : RSNIFTI_OPEN_NONE);
	
    if ( ! p->input->readable ) {
		return;
    }
	
	if ( p->outputpath != NULL ) {
		p->output = fopen(p->outputpath, "w");
	}

	rsSetThreadsNum(p->threads);
	
    // load mask
    if ( p->maskpath != NULL ) {
        p->mask = rsMaskInit(p->maskpath);
        rsMaskLoad(p->mask, p->input);
        if ( ! p->mask->readable) {
            fprintf(stderr, "\nError: Mask (%s) invalid.\n", p->maskpath);
			return;
        }
    }

    if ( p->mask2path != NULL ) {
        p->mask2 = rsMaskInit(p->mask2path);
        rsMaskLoad(p->mask2, p->input);
        if ( ! p->mask2->readable) {
            fprintf(stderr, "\nError: Mask (%s) invalid.\n", p->mask2path);
			return;
        }
    }

    p->parametersValid = TRUE;
}

void rsTimecourseRun(rsTimecourseParameters *p)
{
	p->parametersValid = FALSE;
	
	if ( p->output == NULL ) {
		p->output = stdout;
	}
	
	// extract single voxel timeseries
	if ( p->point != NULL ) {
		rsTimecourseRunSingleVoxelExtraction(p);
		return;
	}

	// remove NaN points
	if ( p->verbose ) {
		fprintf(stdout, "Removing voxels that are either NaN or have a StdDev of 0\n");
	}

	long nRemovedPoints = rsCleanMaskFromNaNs(p->mask, p->input);

	if ( p->verbose ) {
		fprintf(stdout, "Removed %ld voxels, %ld remaining\n", nRemovedPoints, p->mask->nPoints);
	}

	if ( p->mask2 != NULL ) {
		nRemovedPoints = rsCleanMaskFromNaNs(p->mask2, p->input);
		if ( p->verbose ) {
			fprintf(stdout, "Removed %ld voxels, %ld remaining\n", nRemovedPoints, p->mask2->nPoints);
		}
	}

	// run appropriate aggregation algorithm
	switch ( p->algorithm ) {
		case RSTOOLS_TIMECOURSE_ALGORITHM_MEAN:
		case RSTOOLS_TIMECOURSE_ALGORITHM_STDDEV:
			rsTimecourseRunMeanOrStdDev(p);
			break;
		case RSTOOLS_TIMECOURSE_ALGORITHM_SPCA:
		case RSTOOLS_TIMECOURSE_ALGORITHM_TPCA:
			rsTimecourseRunPCA(p);
			break;
		case RSTOOLS_TIMECOURSE_ALGORITHM_CSP:
			rsTimecourseRunCSP(p);
	}
	
	if ( p->outputpath != NULL ) {
		fclose(p->output);
	}
}

void rsTimecourseRunSingleVoxelExtraction(rsTimecourseParameters *p)
{
    // check inputs
    if ( (p->point->x<0) || (p->point->x>=p->input->xDim) ) {
        fprintf(stderr, "\nError: x index (%d) out of range [0..%d]\n", p->point->x, p->input->xDim-1);
		return;
    }
    if ( (p->point->y<0) || (p->point->y>=p->input->yDim) ) {
        fprintf(stderr, "\nError: y index (%d) out of range [0..%d]\n", p->point->y, p->input->yDim-1);
		return;
    }
    if ( (p->point->z<0) || (p->point->z>=p->input->zDim) ) {
        fprintf(stderr, "\nError: z index (%d) out of range [0..%d]\n", p->point->z, p->input->zDim-1);
		return;
    }

	// prepare buffer for one timeseries
    p->input->data = rsMalloc(rsGetBufferSize(1, 1, 1, p->input->vDim, p->input->dt));

    // read out timecourse
    FslReadTimeSeries(p->input->fslio, p->input->data, p->point->x, p->point->y, p->point->z, p->input->vDim);

	// convert and scale
	double series[p->input->vDim];
	convertBufferToScaledDouble(&series[0], p->input->data, p->input->vDim, p->input->slope, p->input->inter, p->input->dt);

	// print
	for (int t=0; t<p->input->vDim; t++) {
		fprintf(p->output, "%.10f\n", series[t]);
	}
	
	p->parametersValid = TRUE;
}

void rsTimecourseRunMeanOrStdDev(rsTimecourseParameters *p)
{
	int t;
	double timecourse[p->input->vDim];

	if ( p->verbose ) {
		if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_MEAN ) {
			fprintf(stdout, "Coputing mean mask timecourse\n");
		} else if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_STDDEV ) {
			fprintf(stdout, "Coputing the standard deviation within the mask for every volume\n");
		}
	}

	// Iterate over all timepoints
	#pragma omp parallel num_threads(rsGetThreadsNum()) private(t) shared(timecourse)
	{
		#pragma omp for schedule(guided)
		for ( t = 0; t<p->input->vDim; t=t+1 ) {

			double sum = 0;

			// Read out datapoints from the buffer
			double *pointValues = (double*)rsMalloc(p->mask->nPoints*sizeof(double));
			rsExtractPointsFromRSNiftiFileBuffer(p->input, pointValues, p->mask->maskPoints, p->mask->nPoints, t);

			// Iterate over all points in the mask and compute mean
			for ( unsigned long i=0; i<p->mask->nPoints; i=i+1) {
				sum = sum + pointValues[i];
			}
			timecourse[t] = sum / p->mask->nPoints;

			// If necessary iterate over them once more to compute the standard deviation
			if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_STDDEV ) {

				sum = 0.0;

				for ( unsigned long i=0; i<p->mask->nPoints; i=i+1) {
					sum = sum + pow(pointValues[i]-timecourse[t], 2.0) / (p->mask->nPoints-1.0);
				}

				timecourse[t] = sqrt(sum);
			}

			free(pointValues);
		}
	}

	// Write out result
	for ( t = 0; t<p->input->vDim; t=t+1 ) {
		fprintf(p->output, "%.10f\n", timecourse[t]);
	}
	
	p->parametersValid = TRUE;
}

void rsTimecourseRunPCA(rsTimecourseParameters *p)
{
	if ( p->verbose ) {
		fprintf(stdout, "Computing PCA components for timecourses in the mask\n");
	}

	// Prepare data matrix
	gsl_matrix* data = gsl_matrix_alloc(p->mask->nPoints, p->input->vDim);
	for ( long n=0; n < p->mask->nPoints; n=n+1) {

		double *timecourse = rsMalloc(sizeof(double) * p->input->vDim);
		rsExtractTimecourseFromRSNiftiFileBuffer(p->input, timecourse, &p->mask->maskPoints[n]);

		double mean  = 0;
		double stdev = 1;

		if ( p->useStandardScores ) {
			mean  = gsl_stats_mean(timecourse, 1, p->input->vDim);
			stdev = gsl_stats_sd(timecourse, 1, p->input->vDim);
		}

		for ( int t = 0; t < p->input->vDim; t=t+1 ) {
			gsl_matrix_set(data, n, t, (timecourse[t] - mean) / stdev);
		}

		free(timecourse);
	}

	// Run PCA
	gsl_matrix* components;
	rsPCAResult *pcaResult;
	long nEigenvalues = p->mask->nPoints;

	if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_SPCA ) {
		pcaResult = rsPCA(data, p->minVariance, p->nComponents, p->verbose);
		components = pcaResult->transformed;

		if ( p->spatialmappath != NULL) {
			rsWriteSpatialMap(p->spatialmappath, p->input, &p->mask->maskPoints[0], pcaResult->eigenvectors);
		}
	} else { // TPCA
		nEigenvalues = p->input->vDim;
		pcaResult = rsTPCA(data, p->minVariance, p->nComponents, p->verbose);
		components = pcaResult->eigenvectors;

		if ( p->spatialmappath != NULL) {
			rsWriteSpatialMap(p->spatialmappath, p->input, &p->mask->maskPoints[0], pcaResult->transformed);
		}
	}
	gsl_matrix_free(data);

	// Write out eigenvalues
	if ( p->eigenvaluespath != NULL ) {
		FILE *fp;
		fp = fopen( p->eigenvaluespath , "w" );

		for ( long i=0; i<nEigenvalues; i=i+1 ) {
			fprintf(fp, "%.10f\n", gsl_vector_get(pcaResult->eigenvalues_all, i));
		}

		fclose(fp);
	}

	if ( p->verbose ) {
		if ( p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_SPCA ) {
			fprintf(stdout, "Reduced data matrix(%dx%d) after spatial PCA:\n", (int)components->size1, (int)components->size2);
		} else {
			fprintf(stdout, "Selected tPCA eigenvectors(%dx%d)\n", (int)components->size1, (int)components->size2);
		}
	}

	for ( int t = 0; t < p->input->vDim; t=t+1 ) {
		for ( int c = 0; c<components->size1; c=c+1) {
           	fprintf(p->output, "% 5.10f", gsl_matrix_get(components, c, t));
			if ( c < (components->size1-1) ) {
				fprintf(p->output, "\t");
			}
		}
		fprintf(p->output, "\n");
	}

	rsPCAResultFree(pcaResult);
	
	p->parametersValid = TRUE;
}

void rsTimecourseRunCSP(rsTimecourseParameters *p)
{
	if ( p->verbose ) {
    	fprintf(stdout, "Computing common spatial patterns for the timecourses in both masks\n");
    }

	if ( p->nComponents < 2 ) {
		p->nComponents = (int) (floor(p->input->vDim/2.0) * 2.0);
	}

    // Prepare data matrices
	long n;
    gsl_matrix* dataA = gsl_matrix_alloc(p->input->vDim, p->mask->nPoints);
    gsl_matrix* dataB = gsl_matrix_alloc(p->input->vDim, p->mask2->nPoints);

    #pragma omp parallel num_threads(rsGetThreadsNum()) shared(p) private(n)
    {
		// first mask
    	#pragma omp for schedule(guided)
    	for ( n=0; n < p->mask->nPoints; n=n+1) {

    		double *timecourse = rsMalloc(sizeof(double) * p->input->vDim);
			rsExtractTimecourseFromRSNiftiFileBuffer(p->input, timecourse, &p->mask->maskPoints[n]);

    		double mean  = 0;
    		double stdev = 1;

    		if ( p->useStandardScores ) {
    			mean  = gsl_stats_mean(timecourse, 1, p->input->vDim);
    			stdev = gsl_stats_sd(timecourse, 1, p->input->vDim);
    		}

    		for ( int t = 0; t < p->input->vDim; t=t+1 ) {
    			gsl_matrix_set(dataA, t, n, (timecourse[t] - mean) / stdev);
    		}

    		free(timecourse);
    	}

		// second mask
    	#pragma omp for schedule(guided)
    	for ( n=0; n < p->mask2->nPoints; n=n+1) {

    		double *timecourse = rsMalloc(sizeof(double) * p->input->vDim);
			rsExtractTimecourseFromRSNiftiFileBuffer(p->input, timecourse, &p->mask2->maskPoints[n]);

    		double mean  = 0;
    		double stdev = 1;

    		if ( p->useStandardScores ) {
    			mean  = gsl_stats_mean(timecourse, 1, p->input->vDim);
    			stdev = gsl_stats_sd(timecourse, 1, p->input->vDim);
    		}

    		for ( int t = 0; t < p->input->vDim; t=t+1 ) {
    			gsl_matrix_set(dataB, t, n, (timecourse[t] - mean) / stdev);
    		}

    		free(timecourse);
    	}
    }

    // Run CSP
    rsCSPResult *cspResult = rsCSP(dataA, dataB, p->nComponents/2, p->verbose);
    long nEigenvalues = cspResult->eigenvalues_all->size;

    if ( p->spatialmappath != NULL) {
    	//rsWriteSpatialMap(spatialmappath, fslio, nonNanPoints, cspResult->eigenvectors);
    }

    gsl_matrix_free(dataA);
    gsl_matrix_free(dataB);

    // Write out eigenvalues
    if ( p->eigenvaluespath != NULL ) {
    	FILE *fp;
    	fp = fopen(p->eigenvaluespath , "w");

    	for ( long i=0; i<nEigenvalues; i=i+1 ) {
    		fprintf(fp, "%.10f\n", gsl_vector_get(cspResult->eigenvalues_all, i));
    	}

    	fclose(fp);
    }

    if ( p->verbose ) {
    	fprintf(stdout, "Selected eigenvectors(%dx%d) of the CSP anaylsis:\n", (int)cspResult->eigenvectors->size1, (int)cspResult->eigenvectors->size2);
    }

	// Print out eigenvectors
    for ( int t = 0; t < p->input->vDim; t=t+1 ) {
    	for ( int c = 0; c < cspResult->eigenvectors->size1; c=c+1) {
        	fprintf(p->output, "% 5.10f", gsl_matrix_get(cspResult->eigenvectors, c, t));
    		if ( c < (cspResult->eigenvectors->size1-1) ) {
    			fprintf(p->output, "\t");
    		}
    	}
    	fprintf(p->output, "\n");
    }

    rsCSPResultFree(cspResult);

	p->parametersValid = TRUE;
}

void rsTimecourseDestroy(rsTimecourseParameters* p)
{
  	if ( p->input != NULL ) {
		rsCloseNiftiFileAndFree(p->input);
		p->input = NULL;
	}

	if ( p->mask != NULL ) {
		rsMaskFree(p->mask);
		p->mask = NULL;
	}

	if ( p->mask2 != NULL ) {
		rsMaskFree(p->mask2);
		p->mask2 = NULL;
	}

	rsTimecourseFreeParams(p);
}

void rsWriteSpatialMap(char *file, const rsNiftiFile *reference, Point3D *points, gsl_matrix *maps)
{
	int x=-1, y=-1, z=-1, t=0;
	short nMaps  = maps->size1;
	long nPoints = maps->size2;

	rsNiftiFile *outputFile = rsCloneNiftiFile(file, reference, RSNIFTI_OPEN_NONE, nMaps);

    if ( ! outputFile->readable ) {
		exit(EXIT_FAILURE);
    }
    FslSetDataType(outputFile->fslio, NIFTI_TYPE_FLOAT32);
	outputFile->dt = NIFTI_TYPE_FLOAT32;
	rsWriteNiftiHeader(outputFile->fslio, "");

	/* prepare buffer */
    outputFile->data = rsMalloc(rsGetBufferSize(outputFile->xDim, outputFile->yDim, outputFile->zDim, outputFile->vDim, NIFTI_TYPE_FLOAT32));
	rsResetRSNiftiFileBufferToValue(outputFile, sqrt(-1.0));

	/* write spatial maps to buffer */
	for (unsigned long p = 0L; p<nPoints; p=p+1L) {
		double data[nMaps];

		for ( short i=0; i<nMaps; i=i+1 ) {
			data[i] = gsl_matrix_get(maps, i, p);
		}
        rsWriteTimecourseToRSNiftiFileBuffer(outputFile, &data[0], &points[p])
    }

	/* write to file */
    FslWriteVolumes(outputFile->fslio, outputFile->data, nMaps);
	rsCloseNiftiFileAndFree(outputFile);
}
