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
        FslSetDataType(fslioResiduals, p.pixtype);
        FslWriteHeader(fslioResiduals);
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
        FslSetDataType(fslioBetas, p.pixtype);
        FslWriteHeader(fslioBetas);
        
        // prepare buffer
        buffsize = (size_t)((size_t)p.xDim*(size_t)p.yDim*(size_t)p.zDim*(size_t)(p.nRegressors+1)*(size_t)p.dt/(size_t)8);
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
        FslSetDataType(fslioFitted, p.pixtype);
        FslWriteHeader(fslioFitted);
        
        // prepare buffer
        buffsize = (size_t)((size_t)p.vDim*(size_t)p.dt/(size_t)8);
        fittedBuffer = malloc(buffsize);
    }
    
    // Prepare empty timecourse
    int emptyValuesLength = p.vDim > p.nRegressors ? p.vDim : p.nRegressors+1;
    double emptybuffer[emptyValuesLength];
    
    for (int i=0; i<emptyValuesLength; i=i+1){
        emptybuffer[i] = log(-1.0);
    }
    
    // Prepare buffer
    buffsize = (size_t)p.xDim*(size_t)p.yDim*(size_t)p.zDim*(size_t)p.vDim*(size_t)p.dt/(size_t)8;
    buffer   = malloc(buffsize);
    
    if (buffer == NULL) {
        fprintf(stdout, "Not enough free memory :-(\n");
        return 1;
    }
    
    FslReadVolumes(p.fslio, buffer, p.vDim);
    short x,y,z;
    double *timecourse;
    double *residuals;
    double *betas;
    double *fitted;
    Point3D point;
    
    #pragma omp parallel num_threads(p.threads) private(y,x,timecourse,residuals,betas,fitted,point) shared(p,emptybuffer,fslioResiduals,buffer,fslioBetas,fslioFitted)
    {
        #pragma omp for schedule(guided)
        for (z=0; z<p.zDim; z=z+1) {
            fprintf(stdout, "Regressing slice Z%03hd/%03hd\n", z+1, p.zDim);        
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
                            // tbd
                        }
                    
                        /* set the fitted value to 0 so that the nifti isn't empty */
                        if ( fslioFitted != NULL ) {
                            // tbd
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
                        p.verbose
                    );
                    
                    rsWriteTimecourseToBuffer(p.fslio, residuals, buffer, p.slope, p.inter, MakePoint3D(x, y, z), p.xDim, p.yDim, p.zDim, p.vDim);
                    
                    free(timecourse);
                    free(residuals);
                    free(betas);
                    free(fitted);
                }
            }
        }
    }
    
    /* Iterate over all voxels that are to be regressed */
//    for (short z=0; z<p.zDim; z=z+1) {
//        fprintf(stdout, "Regressing slice Z%03hd/%03hd\n", z+1, p.zDim);
//        for (short y=0; y<p.yDim; y=y+1) {
//            for (short x=0; x<p.xDim; x=x+1) {
//                
//                /* If it's not in the mask skip it to improve the performance */
//                if (p.mask != NULL && p.mask[z][y][x] < 0.1) {
//                    
//                    /* set the value in the residuals to 0 so that the nifti isn't empty */
//                    if ( fslioResiduals != NULL ) {
//                        rsWriteTimeSeries(fslioResiduals, emptybuffer, x, y, z, p.vDim);
//                    }
//                    
//                    /* set the value in the betas to 0 so that the nifti isn't empty */
//                    if ( fslioBetas != NULL ) {
//                        rsWriteTimeSeries(fslioBetas, emptybuffer, x, y, z, p.nRegressors+1L);
//                    }
//                    
//                    /* set the fitted value to 0 so that the nifti isn't empty */
//                    if ( fslioFitted != NULL ) {
//                        rsWriteTimeSeries(fslioFitted, emptybuffer, x, y, z, p.vDim);
//                    }
//                    
//                    continue;
//                }
//                if (p.verbose) fprintf(stdout, "\n");
//                
//                /* read out timecourse */
//                FslReadTimeSeries(p.fslio, buffer, x, y, z, p.vDim);
//                convertBufferToScaledDouble(signal, buffer, (long)p.vDim, p.slope, p.inter, p.fslio->niftiptr->datatype);
//                
//                /* run the regression */
//                rsLinearRegression(
//                    (int)p.vDim,
//                    signal,
//                    p.nAllRegressors,
//                    (const double**)p.allRegressors,
//                    betas,
//                    residuals,
//                    fitted,
//                    p.verbose
//                );
//                
//                /* write out residuals if desired */
//                if ( p.saveResidualsPath != NULL ) {
//                    convertScaledDoubleToBuffer(fslioResiduals->niftiptr->datatype, residualsBuffer, residuals, p.slope, p.inter, p.vDim, 1, 1, FALSE);
//                    rsWriteTimeSeries(fslioResiduals, residualsBuffer, x, y, z, p.vDim);
//                }
//                
//                /* write out betas if desired */
//                if ( p.saveBetasPath != NULL ) {
//                    convertScaledDoubleToBuffer(fslioBetas->niftiptr->datatype, betasBuffer, betas, p.slope, p.inter, p.nRegressors+1L, 1, 1, FALSE);
//                    rsWriteTimeSeries(fslioBetas, betasBuffer, x, y, z, p.nRegressors+1L);
//                }
//                
//                /* write out fitted values if desired */
//                if ( p.saveFittedPath != NULL ) {
//                    convertScaledDoubleToBuffer(fslioFitted->niftiptr->datatype, fittedBuffer, fitted, p.slope, p.inter, p.vDim, 1, 1, FALSE);
//                    rsWriteTimeSeries(fslioFitted, fittedBuffer, x, y, z, p.vDim);
//                }
//            }
//        }
//    }
    
    if ( p.saveResidualsPath != NULL ) {
        fprintf(stdout, "Write out residuals to: %s\n", p.saveResidualsPath);
        FslWriteVolumes(fslioResiduals, buffer, p.vDim);
        FslClose(fslioResiduals);
        free(fslioResiduals);
    }
    
    if ( p.saveBetasPath != NULL ) {
        FslClose(fslioBetas);
        free(betasBuffer);
    }
    
    if ( p.saveFittedPath != NULL ) {
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

