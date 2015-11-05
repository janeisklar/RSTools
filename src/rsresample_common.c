#include "rsresample_common.h"
#include "rsresample_ui.h"
#include "utils/rsio.h"
#include "maths/utils.h"
#include <math.h>

void rsResampleInit(rsResampleParameters *p)
{
    p->parametersValid = FALSE;
    
    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->outputpath,
        RSIO_LASTFILE
    });

    if ( ! inputsReadable || ! outputsWritable ) {
        return;
    }
    
    rsSetThreadsNum(p->threads);

    if ( p->verbose ) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "Output file: %s\n", p->outputpath);
        fprintf(stdout, "TR(in): %.4f\n", p->inputTR);
        fprintf(stdout, "TR(out): %.4f\n", p->outputTR);
        fprintf(stdout, "Lanczos kernel order: %d\n", p->order);
    }
    
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }
       
    if ( p->verbose ) {
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
    }

    // Determine how many volumes will remain
    int newVDim = (int)floor(((float)p->input->vDim)*(p->inputTR/p->outputTR)-p->order);
        
    if ( p->verbose ) {
        fprintf(stdout, "# of volumes(in): %d\n", p->input->vDim);
        fprintf(stdout, "# of volumes(out): %d\n", newVDim);
    }

    // Create output volume
    p->output = rsCloneNiftiFile(p->outputpath, p->input, RSNIFTI_OPEN_ALLOC, newVDim);
    
    if ( ! p->output->readable ) {
        fprintf(stderr, "\nError: The nifti file containing the resampled output (%s) could not be created.\n", p->outputpath);
        return;
    }
    
    // Check if additional regressors are readable/writable
    for ( int i=0; i<p->nRegressors; i++ ) {
        if ( ! rsCheckInputs((const char*[]){ (const char*)p->regressorInputs[i], RSIO_LASTFILE }) ) {
            return;
        }
        
        if ( ! rsCheckOutputs((const char*[]){ (const char*)p->regressorOutputs[i], RSIO_LASTFILE }) ) {
            return;
        }
        
        if ( p->verbose ) {
            fprintf(stdout, "Addtional regressor file #%d:\n  In: %s\n  Out: %s\n", (i+1), p->regressorInputs[i], p->regressorOutputs[i]);
        }
    }
    
    p->parametersValid = TRUE;
    return;
}

void rsResampleRun(rsResampleParameters *p)
{
    p->parametersValid = FALSE;

    // Run resampling
    const double scaling = p->inputTR/p->outputTR;
    short x,y,z, processedSlices = 0;
    double *signalIn;
    double *signalOut;
    Point3D *point;
    int i,j,k,nResampledValues;
    double **regressors;
    double **resampledRegressor;
    long nRegressors, nValues;
    FILE *fRegressor;
    
    omp_lock_t updateProgressLock;
    omp_init_lock(&updateProgressLock);
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(z,y,x,signalIn,signalOut,point,i,j,k,regressors,nRegressors,nValues,fRegressor,nResampledValues,resampledRegressor) shared(processedSlices)
    {
        // resample additional regressors
        #pragma omp for schedule(guided) nowait
        for (i=0; i < p->nRegressors; i++) {
    
            // load regressor
            regressors = rsLoadRegressors(p->regressorInputs[i], &nRegressors, &nValues, 1.0);
            nResampledValues = (int)floor(((float)nValues)*scaling-p->order);
            resampledRegressor = d2matrix(nRegressors-1, nResampledValues-1);
    
            // resample regressor
            for (j=1; j<=nRegressors; j++) {
                rsInterpolationLanczosConvolve(resampledRegressor[j-1], regressors[j], nValues, nResampledValues, p->order, scaling);
            }
            
            // write to file
            fRegressor = fopen(p->regressorOutputs[i], "w");

            if (fRegressor == NULL) {
                fprintf(stderr, "Error: Regressor file '%s' could not be saved (wrong permissions?).\n", p->regressorOutputs[i]);
            }
            

            for (k=0; k < nResampledValues; k++) {
                for (j=0; j < nRegressors; j++) {
                    if ( j > 0 ) {
                        fprintf(fRegressor, " ");
                    }
                    fprintf(fRegressor, "%.10f", resampledRegressor[j][k]);
                }
                if ( k+1 < nResampledValues ) {
                    fprintf(fRegressor, "\n");
                }
            }
            
            fclose(fRegressor);
            
            // free space
            free(resampledRegressor[0]);
            free(resampledRegressor);
        }

        // resample nifti
        #pragma omp for schedule(guided) nowait
        for (z=0; z<p->input->zDim; z++) {
            for (y=0; y<p->input->yDim; y++) {
                for (x=0; x<p->input->xDim; x++) {
                    
                    point = rsMakePoint3D(x, y, z);
                    
                    // read out timecourse from buffer
                    signalIn  = (double*)rsMalloc(p->input->vDim*sizeof(double));
                    
                    rsExtractTimecourseFromRSNiftiFileBuffer(p->input, signalIn, point);
                    
                    // apply filter
                    signalOut = (double*)rsMalloc(p->output->vDim*sizeof(double));
                    rsInterpolationLanczosConvolve(signalOut, signalIn, p->input->vDim, p->output->vDim, p->order, scaling);
                    
                    // write out resampled data to buffer
                    rsWriteTimecourseToRSNiftiFileBuffer(p->output, signalOut, point);

                    free(signalIn);
                    free(signalOut);
                }
            }
            
            /* show progress */
            if (p->progressCallback != NULL) {
                omp_set_lock(&updateProgressLock);
                rsReportProgressEvent *event = (rsReportProgressEvent*)rsMalloc(sizeof(rsReportProgressEvent));
                event->run = processedSlices;
                processedSlices += 1;
                event->percentage = (double)processedSlices*100.0 / (double)p->input->zDim;
                rsReportProgressCallback_t cb = p->progressCallback->cb;
                void *data = p->progressCallback->data;
                cb(event, data);
                rsFree(event);
                omp_unset_lock(&updateProgressLock);
            } else if (p->verbose) {
                #pragma omp atomic
                processedSlices += 1;
            
                if (processedSlices > 0 && processedSlices % (short)(p->input->vDim >= 10 ? p->input->vDim / 10 : p->input->vDim) == 0) {
                    fprintf(stdout, "..%.0f%%\n", ceil((float)processedSlices*100.0 / (float)p->input->zDim));
                }
            }
        }
    }

    omp_destroy_lock(&updateProgressLock);
    
    if ( p->verbose ) {
        fprintf(stdout, "Write out result to: %s\n", p->outputpath);
    }
    
    rsWriteNiftiHeader(p->output->fslio, p->callString);
    FslWriteVolumes(p->output->fslio, p->output->data, p->output->vDim);

    p->parametersValid = TRUE;
}

void rsResampleDestroy(rsResampleParameters *p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFileAndFree(p->input);
        p->input = NULL;
    }
    
    if ( p->output != NULL ) {
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }
    
    rsResampleFreeParams(p);
}
