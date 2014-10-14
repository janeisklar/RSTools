#include "rsregression_common.h"
#include "utils/rsio.h"

void rsRegressionInit(rsRegressionParameters* p)
{
    p->parametersValid = FALSE;

    /* verify accessibility of inputs/outüuts */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        (const char*)p->regressorspath,
        (const char*)p->maskpath,
        (const char*)p->regressorCommentPath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->saveResidualsPath,
        (const char*)p->saveFittedPath,
        (const char*)p->saveBetasPath,
        (const char*)p->savemaskpath,
        RSIO_LASTFILE
    });
    
    if ( ! inputsReadable || ! outputsWritable ) {
        return;
    }    
    
    rsSetThreadsNum(p->threads);
    
    p->filterActive = p->TR >= 0.0 || p->freqLow >= 0.0 || p->freqHigh >= 0.0;
    
    if ( p->verbose ) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "Regressors file: %s\n", p->regressorspath);
        fprintf(stdout, "Mask file: %s\n", p->maskpath);
        fprintf(stdout, "Residuals file: %s\n", p->saveResidualsPath);
        fprintf(stdout, "Fitted file: %s\n", p->saveFittedPath);
        fprintf(stdout, "Betas file: %s\n", p->saveBetasPath);
        
        if ( p->filterActive ) {
            fprintf(stdout, "Bandpass filter active!\n");
            fprintf(stdout, "Sampling rate: %.4f\n", p->TR);
            fprintf(stdout, "Bandpass: (%.4f-%.4fHz)\n", p->freqLow, p->freqHigh);
        }
    }
    
    // load regressors
    p->regressors = rsLoadRegressors(p->regressorspath, &p->nRegressors, &p->nRegressorValues, 1.0);

    if ( p->verbose ) {
        fprintf(stdout, "Regressors: %ld, Samples: %ld\n", p->nRegressors, p->nRegressorValues);
    }
    
    // open input file
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }

    // create residuals file (output)
    if ( p->saveResidualsPath != NULL ) {
        
        p->residuals = rsCloneNiftiFile(p->saveResidualsPath, p->input, RSNIFTI_CLONE_POINTER, RSNIFTI_CLONE_AS_INPUT);

        if ( ! p->residuals->readable ) {
            fprintf(stderr, "\nError: The nifti file containing the residuals output (%s) could not be created.\n", p->saveResidualsPath);
            return;
        }       
    }

    // create betas file (output)
    if ( p->saveBetasPath != NULL ) {
        
        p->betas = rsCloneNiftiFile(p->saveBetasPath, p->input, RSNIFTI_OPEN_ALLOC, p->nRegressors+1);

        if ( ! p->betas->readable ) {
            fprintf(stderr, "\nError: The nifti file containing the betas output (%s) could not be created.\n", p->saveBetasPath);
            return;
        }       
    }
    
    // create fitted values file (output)
    if ( p->saveFittedPath != NULL ) {
        
        p->fitted = rsCloneNiftiFile(p->saveFittedPath, p->input, RSNIFTI_OPEN_ALLOC, RSNIFTI_CLONE_AS_INPUT);

        if ( ! p->fitted->readable ) {
            fprintf(stderr, "\nError: The nifti file containing the fitted values (%s) could not be created.\n", p->saveFittedPath);
            return;
        }       
    }
    
    // prepare bandpass filter
    if ( p->filterActive ) {
    
        p->nyquist_frequency = 1.0 / (2.0 * p->TR);
        p->bin_width         = p->nyquist_frequency / (p->input->vDim/2);
        
        // Compute the number of frequency regressors that will be added to the existing regressors
        p->nFrequencyBinsLow    = (int)floor(p->freqLow / p->bin_width);                           // number of frequency bins before the lower cutoff frequency
        p->nFrequencyBinsHigh   = (int)floor((p->nyquist_frequency - p->freqHigh) / p->bin_width); // number of frequency bins after the highter cutoff frequency
        p->nFrequencyBins       = p->nFrequencyBinsLow + p->nFrequencyBinsHigh;                    // number of frequency bins in total
        p->nFrequencyRegressors = p->nFrequencyBins * 2;                                           // number of frequency regressors (both sine and cosine)
        
        // Compute the frequencies for the bins
        
        if (p->verbose) 
            fprintf(stdout, "Bin width: %.4f\n", p->bin_width);
            
        p->frequencyBins = rsMalloc(p->nFrequencyBins * sizeof(double));
        
        for ( int i=0; i<p->nFrequencyBins; i=i+1 ) {
            if ( (i < p->nFrequencyBinsLow) ) {
                p->frequencyBins[i] = (i + 1) * p->bin_width;
            } else {
                p->frequencyBins[i] = (i - p->nFrequencyBinsLow + 1) * p->bin_width + p->freqHigh;
            }
            
            if (p->verbose)
                fprintf(stdout, "Frequency bin %d(%s): %.4f\n", i, (i < p->nFrequencyBinsLow ? "low" : "high"), p->frequencyBins[i]);
        }
    }
    
    // Add frequency regressors to the other regressors if requested
    p->nAllRegressors = p->filterActive ? p->nRegressors + 1 + p->nFrequencyRegressors : p->nRegressors + 1;
    
    if ( p->filterActive ) {
        p->allRegressors = d2matrix(p->nAllRegressors-1, p->input->vDim-1);
        
        for (int i=0; i<p->nAllRegressors; i=i+1) {
            
            if ( i < p->nRegressors ) {
                for (int t=0; t<p->input->vDim; t=t+1) {
                    p->allRegressors[i][t] = p->regressors[i][t];
                }
            } else if ( i < (p->nRegressors + p->nFrequencyBins) ) {
                const int j = i - p->nRegressors;
                for (int t=0; t<p->input->vDim; t=t+1) {
                    p->allRegressors[i][t] = rsSampleSineWave(p->TR, p->frequencyBins[j], t);
                }
            } else {
                const int j = i - p->nRegressors - p->nFrequencyBins;
                for (int t=0; t<p->input->vDim; t=t+1) {
                    p->allRegressors[i][t] = rsSampleCosineWave(p->TR, p->frequencyBins[j], t);
                }
            }
        }
    } else {
        p->allRegressors = p->regressors;
    }

    // if an additional regressor comment file was specified load it
    if ( p->regressorCommentPath != NULL ) {
        p->regressorComment = rsReadCommentFile(p->regressorCommentPath);
    }

    // assemble the comment that will be stored in the output niftis
    char *comment;
    if ( p->regressorComment == NULL ) {
        comment = p->callString;
    } else {
        char *separator = "\nnRegressor Info:\n";
        comment = malloc(sizeof(char)*(strlen(p->callString)+strlen(separator)+strlen(p->regressorComment)+1));
        sprintf(&comment[0], "%s%s%s", p->callString, separator, p->regressorComment);
    }
    p->comment = comment;
    
    // load mask
    p->mask = NULL;
    if ( p->maskpath != NULL ) {
        unsigned long nPoints = 0L;
        p->mask = d3matrix(p->input->zDim-1, p->input->yDim-1, p->input->xDim-1);
        Point3D *maskPoints = rsReadMask(p->maskpath, p->input->xDim, p->input->yDim, p->input->zDim, &nPoints, p->savemaskpath, p->input->fslio, p->mask);
        if ( maskPoints == NULL) {
            fprintf(stderr, "\nError: Mask invalid.\n");
            return;
        }
        free(maskPoints);
    }
    
    p->parametersValid = TRUE;    
}

void rsRegressionRun(rsRegressionParameters *p)
{
    
    // Prepare empty timecourse (containing NaNs)
    int emptyValuesLength = p->input->vDim > p->nRegressors ? p->input->vDim : p->nRegressors+1;
    double emptybuffer[emptyValuesLength];
    
    for (int i=0; i<emptyValuesLength; i=i+1){
        emptybuffer[i] = log(-1.0);
    }
    
    short x,y,z, processedSlices = 0;
    double *timecourse;
    double *residuals;
    double *betas;
    double *fitted;
    Point3D *point;
    omp_lock_t updateProgressLock;
    omp_init_lock(&updateProgressLock);
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(y,x,timecourse,residuals,betas,fitted,point) shared(p,emptybuffer)
    {
        #pragma omp for schedule(guided)
        for (z=0; z<p->input->zDim; z++) {
            for (y=0; y<p->input->yDim; y++) {
                for (x=0; x<p->input->xDim; x++) {
                    
                    point = rsMakePoint3D(x, y, z);
                    
                    // If it's not in the mask skip it to improve the performance
                    if (p->mask != NULL && p->mask[z][y][x] < 0.1) {
                    
                        // set the value in the residuals to NaN so that the nifti isn't empty
                        if ( p->residuals != NULL ) {
                            rsWriteTimecourseToRSNiftiFileBuffer(p->residuals, emptybuffer, point);
                        }
                    
                        // set the value in the betas to NaN so that the nifti isn't empty
                        if ( p->betas != NULL ) {
                            rsWriteTimecourseToRSNiftiFileBuffer(p->betas, emptybuffer, point);
                        }
                    
                        // set the fitted value to NaN so that the nifti isn't empty
                        if ( p->fitted != NULL ) {
                            rsWriteTimecourseToRSNiftiFileBuffer(p->fitted, emptybuffer, point);
                        }
                        
                        free(point);
                        continue;
                    }
                    
                    timecourse = rsMalloc(p->input->vDim*sizeof(double));
                    residuals  = rsMalloc(p->input->vDim*sizeof(double));
                    betas      = rsMalloc(p->nAllRegressors*sizeof(double));
                    fitted     = rsMalloc(p->input->vDim*sizeof(double));
                    
                    rsExtractTimecourseFromRSNiftiFileBuffer(p->input, timecourse, point);
                    
                    // run the regression
                    rsLinearRegression(
                        (int)p->input->vDim,
                        timecourse,
                        p->nAllRegressors,
                        (const double**)p->allRegressors,
                        betas,
                        residuals,
                        fitted,
                        p->zScoreRegression,
                        p->verbose
                    );

                    // write the results to the buffers
                    if ( p->residuals != NULL ) {
                        rsWriteTimecourseToRSNiftiFileBuffer(p->residuals, residuals, point);
                    }

                    if ( p->betas != NULL ) {
                        rsWriteTimecourseToRSNiftiFileBuffer(p->betas, betas, point);
                    }

                    if ( p->fitted != NULL ) {
                        rsWriteTimecourseToRSNiftiFileBuffer(p->fitted, fitted, point);
                    }
                    
                    free(timecourse);
                    free(residuals);
                    free(betas);
                    free(fitted);
                    free(point);
                }
            }
            
            /* show progress */
            if ( p->progressCallback != NULL ) {
                omp_set_lock(&updateProgressLock);
                rsReportProgressEvent *event = (rsReportProgressEvent*)rsMalloc(sizeof(rsReportProgressEvent));
                event->run = processedSlices;
                processedSlices += 1;
                event->percentage = ((double)processedSlices*100.0) / (double)p->input->zDim;
                rsReportProgressCallback_t cb = p->progressCallback->cb;
                void *data = p->progressCallback->data;
                cb(event, data);
                rsFree(event);
                omp_unset_lock(&updateProgressLock);
            } else if ( p->verbose ) {
                #pragma omp atomic
                processedSlices += 1;
            
                if (processedSlices > 0 && processedSlices % (short)(p->input->zDim / 10) == 0) {
                    fprintf(stdout, "..%.0f%%\n", ceil((float)processedSlices*100.0 / (float)p->input->zDim));
                }
            }
        }
    }

    omp_destroy_lock (&updateProgressLock);
    
    /* Write out buffers to the corresponding files */
    if ( p->residuals != NULL ) {
        if (p->verbose) fprintf(stdout, "Saving residuals.. (%s)\n", p->saveResidualsPath);
        rsWriteNiftiHeader(p->residuals->fslio, p->comment);
        FslWriteVolumes(p->residuals->fslio, p->residuals->data, p->residuals->vDim);
    }
    
    if ( p->betas != NULL ) {
        if (p->verbose) fprintf(stdout, "Saving betas.. (%s)\n", p->saveBetasPath);
        rsWriteNiftiHeader(p->betas->fslio, p->comment);
        FslWriteVolumes(p->betas->fslio, p->betas->data, p->betas->vDim);
    }
    
    if ( p->fitted != NULL ) {
        if (p->verbose) fprintf(stdout, "Saving fitted values.. (%s)\n", p->saveFittedPath);
        rsWriteNiftiHeader(p->fitted->fslio, p->comment);
        FslWriteVolumes(p->fitted->fslio, p->fitted->data, p->fitted->vDim);
    }
}

void rsRegressionDestroy(rsRegressionParameters* p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFileAndFree(p->input);
        p->input = NULL;
    }
    
    if ( p->residuals != NULL ) {
        p->residuals->data = NULL;
        rsCloseNiftiFileAndFree(p->residuals);
        p->residuals = NULL;
    }
    
    if ( p->fitted != NULL ) {
        rsCloseNiftiFileAndFree(p->fitted);
        p->fitted = NULL;
    }

    if ( p->betas != NULL ) {
        rsCloseNiftiFileAndFree(p->betas);
        p->betas = NULL;
    }
    
    if ( p->maskpath != NULL ) {
        free(p->mask);
        p->mask = NULL;
    }

    rsRegressionFreeParams(p);
}
