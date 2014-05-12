/*******************************************************************
 *
 * rsregression2.c
 *
 * rsregression2: Given a 4D-Nifti and a txt file with regressors(columns) this tool will perform a multiple linear regression on it
 *
 * Usage: rsregression2 -input <volume> -regressors <txtFile> [-residuals <volume> | -fitted <volume> | -betas <volume> | -mask <volume>]
 *
 *
 * André Hoffmann
 *******************************************************************/

#include <omp.h>
#include "rsregression_common.h"

int show_help( )
{
    printf(
	   RSTOOLS_VERSION_LABEL "\n"
        "rsregression: Given a 4D-Nifti and a txt file with regressors(columns),\n"
        "              this tool will perform a multiple linear regression on it.\n"
        "\n"
    );
    
    rsRegressionPrintHelp();
    
    return 0;
}

int main(int argc, char * argv[])
{
	if( argc < 2 ) return show_help();
    
    struct rsRegressionParameters p = rsRegressionLoadParams(argc, argv);
    
    if ( ! p.parametersValid ) {
        fprintf(stderr, "Invalid parameters\n");
        return 1;
    }
    
    void *buffer;
	size_t buffsize;
	char *callString = rsMergeStringArray(argc, argv);
	char *comment;
	if ( p.regressorComment == NULL ) {
		comment = callString;
	} else {
		char *separator = "\nRegressor Info:\n";
		comment = malloc(sizeof(char)*(strlen(callString)+strlen(separator)+strlen(p.regressorComment)+1));
		sprintf(&comment[0], "%s%s%s", callString, separator, p.regressorComment);
		free(callString);
	}

    /* prepare residuals file */
    FSLIO *fslioResiduals = NULL;
    
    if ( p.saveResidualsPath != NULL ) {
        
        fslioResiduals = FslOpen(p.saveResidualsPath, "wb");
        
        if (fslioResiduals == NULL) {
            fprintf(stderr, "\nError, could not read header info for %s.\n", p.saveResidualsPath);
            return 1;
        }
        
        FslCloneHeader(fslioResiduals, p.fslio);
        FslSetDim(fslioResiduals, p.xDim, p.yDim, p.zDim, p.vDim);
        FslSetDimensionality(fslioResiduals, 4);
        FslSetDataType(fslioResiduals, p.dt);
		rsWriteNiftiHeader(fslioResiduals, comment);
    }
    
    /* prepare betas file */
    FSLIO *fslioBetas = NULL;
   	void *betasBuffer;
    
    if ( p.saveBetasPath != NULL ) {
        
        fslioBetas = FslOpen(p.saveBetasPath, "wb");
        
        if (fslioBetas == NULL) {
            fprintf(stderr, "\nError, could not read header info for %s.\n", p.saveBetasPath);
            return 1;
        }
        
        FslCloneHeader(fslioBetas, p.fslio);
        FslSetDim(fslioBetas, p.xDim, p.yDim, p.zDim, (p.nRegressors+1L));
        FslSetDimensionality(fslioBetas, 4);
        FslSetDataType(fslioBetas, p.dt);
		rsWriteNiftiHeader(fslioBetas, comment);
        
        // prepare buffer
		buffsize = rsGetBufferSize(p.xDim, p.yDim, p.zDim, p.nRegressors+1, p.dt);
        betasBuffer = malloc(buffsize);
    }
    
    /* prepare fitted values file */
    FSLIO *fslioFitted = NULL;
   	void *fittedBuffer;
    
    if ( p.saveFittedPath != NULL ) {
        
        fslioFitted = FslOpen(p.saveFittedPath, "wb");
        
        if (fslioFitted == NULL) {
            fprintf(stderr, "\nError, could not read header info for %s.\n",p.saveFittedPath);
            return 1;
        }
        
        FslCloneHeader(fslioFitted, p.fslio);
        FslSetDim(fslioFitted, p.xDim, p.yDim, p.zDim, p.vDim);
        FslSetDimensionality(fslioFitted, 4);
        FslSetDataType(fslioFitted, p.dt);
		rsWriteNiftiHeader(fslioFitted, comment);
        
        // prepare buffer
		buffsize = rsGetBufferSize(p.xDim, p.yDim, p.zDim, p.vDim, p.dt);
        fittedBuffer = malloc(buffsize);
    }
	free(callString);
    
    // Prepare empty timecourse
    int emptyValuesLength = p.vDim > p.nRegressors ? p.vDim : p.nRegressors+1;
    double emptybuffer[emptyValuesLength];
    
    for (int i=0; i<emptyValuesLength; i=i+1){
        emptybuffer[i] = log(-1.0);
    }
    
    // Prepare buffer
	buffsize = rsGetBufferSize(p.xDim, p.yDim, p.zDim, p.vDim, p.dt);
    buffer   = malloc(buffsize);
        
    if (buffer == NULL) {
        fprintf(stdout, "Not enough free memory :-(\n");
        return 1;
    }
    
    FslReadVolumes(p.fslio, buffer, p.vDim);
    short x,y,z, processedSlices = 0;
    double *timecourse;
    double *residuals;
    double *betas;
    double *fitted;
    Point3D point;
    
    #pragma omp parallel num_threads(p.threads) private(y,x,timecourse,residuals,betas,fitted,point) shared(p,emptybuffer,fslioResiduals,buffer,fslioBetas,fslioFitted,betasBuffer,fittedBuffer)
    {
        #pragma omp for schedule(guided)
        for (z=0; z<p.zDim; z=z+1) {
            for (y=0; y<p.yDim; y=y+1) {
                for (x=0; x<p.xDim; x=x+1) {
                    
                    point = MakePoint3D(x, y, z);
                    
                    /* If it's not in the mask skip it to improve the performance */
                    if (p.mask != NULL && p.mask[z][y][x] < 0.1) {
                    
                        /* set the value in the residuals to 0 so that the nifti isn't empty */
                        if ( fslioResiduals != NULL ) {
                            rsWriteTimecourseToBuffer(fslioResiduals, emptybuffer, buffer, p.slope, p.inter, point, p.xDim, p.yDim, p.zDim, p.vDim);
                        }
                    
                        /* set the value in the betas to 0 so that the nifti isn't empty */
                        if ( fslioBetas != NULL ) {
                            rsWriteTimecourseToBuffer(fslioBetas, emptybuffer, betasBuffer, p.slope, p.inter, point, p.xDim, p.yDim, p.zDim, p.nRegressors+1L);
                        }
                    
                        /* set the fitted value to 0 so that the nifti isn't empty */
                        if ( fslioFitted != NULL ) {
                            rsWriteTimecourseToBuffer(fslioFitted, emptybuffer, fittedBuffer, p.slope, p.inter, point, p.xDim, p.yDim, p.zDim, p.vDim);
                        }
                        
                        continue;
                    }
                    
                    timecourse = malloc(p.vDim*sizeof(double));
                    residuals  = malloc(p.vDim*sizeof(double));
                    betas      = malloc(p.nAllRegressors*sizeof(double));
                    fitted     = malloc(p.vDim*sizeof(double));
                    
                    rsExtractTimecourseFromBuffer(p.fslio, timecourse, buffer, p.slope, p.inter, MakePoint3D(x, y, z), p.xDim, p.yDim, p.zDim, p.vDim);
                    
                    /* run the regression */
                    rsLinearRegression(
                        (int)p.vDim,
                        timecourse,
                        p.nAllRegressors,
                        (const double**)p.allRegressors,
                        betas,
                        residuals,
                        fitted,
                        p.zScoreRegression,
                        p.verbose
                    );
                    
                    if ( fslioResiduals != NULL ) {
                        rsWriteTimecourseToBuffer(fslioResiduals, residuals, buffer, p.slope, p.inter, MakePoint3D(x, y, z), p.xDim, p.yDim, p.zDim, p.vDim);
                    }
                    
                    if ( fslioBetas != NULL ) {
                        rsWriteTimecourseToBuffer(fslioBetas, betas, betasBuffer,  p.slope, p.inter, MakePoint3D(x, y, z), p.xDim, p.yDim, p.zDim, p.nRegressors+1);
                    }
                    
                    if ( fslioFitted != NULL ) {
                        rsWriteTimecourseToBuffer(fslioFitted, fitted, fittedBuffer, p.slope, p.inter, MakePoint3D(x, y, z), p.xDim, p.yDim, p.zDim, p.vDim);
                    }
                    
                    free(timecourse);
                    free(residuals);
                    free(betas);
                    free(fitted);
                }
            }
            
			/* show progress */
			if (p.verbose) {
            	#pragma omp atomic
            	processedSlices += 1;
            
            	if (processedSlices > 0 && processedSlices % (short)(p.zDim / 10) == 0) {
                	fprintf(stdout, "..%.0f%%\n", ceil((float)processedSlices*100.0 / (float)p.zDim));
            	}
			}
        }
    }
    
    /* Write out buffers to the corresponding files */
    if ( p.saveResidualsPath != NULL ) {
        if (p.verbose) fprintf(stdout, "Write out residuals to: %s\n", p.saveResidualsPath);
        FslWriteVolumes(fslioResiduals, buffer, p.vDim);
        FslClose(fslioResiduals);
        free(fslioResiduals);
    }
    
    if ( p.saveBetasPath != NULL ) {
        if (p.verbose) fprintf(stdout, "Write out betas to: %s\n", p.saveBetasPath);
        FslWriteVolumes(fslioBetas, betasBuffer, p.nRegressors+1);
        FslClose(fslioBetas);
        free(betasBuffer);
    }
    
    if ( p.saveFittedPath != NULL ) {
        if (p.verbose) fprintf(stdout, "Write out fitted values to: %s\n", p.saveFittedPath);
        FslWriteVolumes(fslioFitted, fittedBuffer, p.vDim);
        FslClose(fslioFitted);
        free(fittedBuffer);
    }
    
    if ( p.maskpath != NULL ) {
        free(p.mask);
    }
    
    free(buffer);
    FslClose(p.fslio);
    free(p.regressors[0]);
    free(p.regressors);
    if ( p.filterActive ) {
        free(p.allRegressors[0]);
        free(p.allRegressors);
        free(p.frequencyBins);
    }
        
	return 0;
}

