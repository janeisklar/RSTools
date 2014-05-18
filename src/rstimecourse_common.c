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

    // open input file
	BOOL fileNeedsToBeRead = p->algorithm != RSTOOLS_TIMECOURSE_ALGORITHM_MEAN && 
	                         p->algorithm != RSTOOLS_TIMECOURSE_ALGORITHM_STDDEV;
	
	p->input = rsOpenNiftiFile(p->inputpath, fileNeedsToBeRead ? RSNIFTI_OPEN_READ : RSNIFTI_OPEN_NONE);
    if ( ! p->input->readable ) {
		return;
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
	// extract single voxel timeseries
	if ( p->point != NULL ) {
		rsTimecourseRunSingleVoxelExtraction(p);
		return;
	}
	
	// remove NaN points
	if ( p->verbose ) fprintf(stdout, "Removing voxels that are either NaN or have a StdDev of 0\n");
	long nRemovedPoints = rsCleanMaskFromNaNs(p->mask, p->input);
	if ( p->verbose ) fprintf(stdout, "Removed %ld voxels, %ld remaining\n", nRemovedPoints, p->mask->nPoints);
	
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
		if ( p->verbose ) {
            fprintf(stdout, "%d: %.10f\n", t, series[t]);
        } else {
            fprintf(stdout, "%.10f\n", series[t]);
        }
	}
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
		fprintf(stdout, "%.10f\n", timecourse[t]);
	}
}

void rsTimecourseRunPCA(rsTimecourseParameters *p)
{
	
}

void rsTimecourseRunCSP(rsTimecourseParameters *p)
{
	
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
	//rsResetBufferToValue(outputFile->dt, outputFile->data, outputFile->slope, outputFile->inter, outputFile->xDim, outputFile->yDim, outputFile->zDim, nMaps, sqrt(-1.0));
	
	/* write spatial maps to buffer */
	for (unsigned long p = 0L; p<nPoints; p=p+1L) {
		double data[nMaps];
		
		for ( short i=0; i<nMaps; i=i+1 ) {
			data[i] = gsl_matrix_get(maps, i, p);
		}
        rsWriteTimecourseToRSNiftiFileBuffer(outputFile, &data[0], &points[p])
        //rsWriteTimecourseToBuffer(outputFile->dt, &data[0], outputFile->data, outputFile->slope, outputFile->inter, &points[p], outputFile->xDim, outputFile->yDim, outputFile->zDim, nMaps);
    }
    
	/* write to file */
    FslWriteVolumes(outputFile->fslio, outputFile->data, nMaps);
	rsCloseNiftiFileAndFree(outputFile);
}