#include "rsbandpass_common.h"
#include "rsbandpass_ui.h"
#include "utils/rsio.h"

void rsBandpassInit(rsBandpassParameters *p)
{
    p->parametersValid = FALSE;
    
    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        (const char*)p->maskpath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->saveFilteredPath,
        (const char*)p->savemaskpath,
        (const char*)p->saveAttenuationPath,
        RSIO_LASTFILE
    });

    if ( ! inputsReadable || ! outputsWritable ) {
        return;
    }
    
    rsSetThreadsNum(p->threads);

    if ( p->verbose ) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "Mask file: %s\n", p->maskpath);
        fprintf(stdout, "Filtered file: %s\n", p->saveFilteredPath);
        fprintf(stdout, "lower freq.: %.4f\n", p->freqLow);
        fprintf(stdout, "upper freq.: %.4f\n", p->freqHigh);
        fprintf(stdout, "TR: %.4f\n", p->TR);
    }

    rsFFTSetEngine(RSFFTFILTER_ENGINE_GSL);
#if RS_FFTW_ENABLED == 1
    if ( p->fftw ) rsFFTSetEngine(RSFFTFILTER_ENGINE_FFTW);
#endif
    
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }
       
    if ( p->verbose ) {
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
    }

    if ( p->paddedT == 0 ) {
        p->paddedT = p->input->vDim;
    } else if ( p->verbose ) {
        fprintf(stdout, "Padding data to have a sampling length of %ld\n", p->paddedT);
    }
    
    if ( p->input->vDim > p->paddedT ) {
        fprintf(stderr, "\nError: datalength(%ld) needs to be longer or equal to the temporal length(%d) of the supplied nifti file: %s.\n",p->paddedT, p->input->vDim, p->inputpath);
        return;
    }
    
    // Load mask
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
    
    // Prepare FFT filter
    p->fftParams = rsFFTFilterInit(p->input->vDim, p->paddedT, p->TR, p->freqLow, p->freqHigh, p->rolloff_method, p->rolloff, p->keepMean, p->verbose);
    
    if ( p->saveAttenuationPath != NULL ) {
        if ( p->verbose ) {
            fprintf(stdout, "Writing attenuation weights to: %s\n", p->saveAttenuationPath);
        }
        FILE *file;
        file = fopen(p->saveAttenuationPath, "wb");
        
        for ( int i=0; i<p->fftParams->paddedT; i=i+1 ) {
            fprintf(file,"%.10f\t%.10f\n", p->fftParams->frequencyBins[i], p->fftParams->binAttenuation[i]);
        }
        
        fclose(file);
    }

    // Create output volume
    p->filteredOutput = rsCloneNiftiFile(p->saveFilteredPath, p->input, RSNIFTI_CLONE_POINTER, RSNIFTI_CLONE_AS_INPUT);
    
    if ( ! p->filteredOutput->readable ) {
        fprintf(stderr, "\nError: The nifti file containing the filtered output (%s) could not be created.\n", p->saveFilteredPath);
        return;
    }
    
    p->parametersValid = TRUE;
    return;
}

void rsBandpassRun(rsBandpassParameters *p)
{
    p->parametersValid = FALSE;
    
    // Prepare empty timecourse
    double emptybuffer[p->input->vDim];
    
    for (int i=0; i<p->input->vDim; i=i+1){
        emptybuffer[i] = log(-1.0);
    }

    // Run FFT bandpass filter
    short x,y,z, processedSlices = 0;
    double *signal;
    Point3D *point;
    
    omp_lock_t updateProgressLock;
    omp_init_lock(&updateProgressLock);
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(z,y,x,signal,point) shared(emptybuffer,processedSlices)
    {
        #pragma omp for schedule(guided)
        for (z=0; z<p->input->zDim; z++) {
            for (y=0; y<p->input->yDim; y++) {
                for (x=0; x<p->input->xDim; x++) {
                    
                    point = rsMakePoint3D(x, y, z);
                    
                    /* If it's not in the mask skip it to improve the performance */
                    if (p->mask != NULL && p->mask[z][y][x] < 0.1) {
                    
                        /* set the value in the filtered data to NaN so that the nifti isn't empty */
                        rsWriteTimecourseToRSNiftiFileBuffer(p->input, emptybuffer, point);
                        continue;
                    }
                    
                    /* read out timecourse from buffer */
                    signal = rsMalloc(p->input->vDim*sizeof(double));
                    rsExtractTimecourseFromRSNiftiFileBuffer(p->input, signal, point);
                    
                    /* apply filter */
                    rsFFTFilter(p->fftParams, signal);
                    
                    /* write out filtered data to buffer */
                    rsWriteTimecourseToRSNiftiFileBuffer(p->filteredOutput, signal, point);

                    free(signal);
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
            
                if (processedSlices > 0 && processedSlices % (short)(p->input->zDim / 10) == 0) {
                    fprintf(stdout, "..%.0f%%\n", ceil((float)processedSlices*100.0 / (float)p->input->zDim));
                }
            }
        }
    }

    omp_destroy_lock(&updateProgressLock);
    
    if ( p->verbose ) {
        fprintf(stdout, "Write out result to: %s\n", p->saveFilteredPath);
    }
    
    rsWriteNiftiHeader(p->filteredOutput->fslio, p->callString);
    FslWriteVolumes(p->filteredOutput->fslio, p->filteredOutput->data, p->filteredOutput->vDim);

    p->parametersValid = TRUE;
}

void rsBandpassDestroy(rsBandpassParameters *p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFileAndFree(p->input);
        p->input = NULL;
    }
    
    if ( p->filteredOutput != NULL ) {
        p->filteredOutput->data = NULL;
        rsCloseNiftiFileAndFree(p->filteredOutput);
        p->filteredOutput = NULL;
    }
    
    if ( p->maskpath != NULL ) {
        free(p->mask);
        p->mask = NULL;
    }

    rsFFTFilterFree(p->fftParams);
    p->fftParams = NULL;
    
    rsBandpassFreeParams(p);
}
