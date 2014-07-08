#include "rscorrelation_common.h"

void rsCorrelationInit(rsCorrelationParameters* p)
{
    p->parametersValid = FALSE;

    /* check if the required arguments have been provided */
    if ( p->inputpath == NULL ) {
        fprintf(stderr, "No input volume specified!\n");
        return;
    }
    
    if ( p->outputpath == NULL ) {
        fprintf(stderr, "No output volume specified!\n");
        return;
    }
    
    rsSetThreadsNum(p->threads);
    
    /* if an additional comment file was specified load it */
    if ( p->commentpath != NULL ) {
        p->comment = rsReadCommentFile(p->commentpath);
    }
    
    /* assemble the comment that will be stored in the resulting output nifti */
    char *comment;
    if ( p->comment == NULL ) {
        comment = p->callString;
    } else {
        char *separator = "\nReference Timecourse Info:\n";
        comment = malloc(sizeof(char)*(strlen(p->callString)+strlen(separator)+strlen(p->comment)+1));
        sprintf(&comment[0], "%s%s%s", p->callString, separator, p->comment);
    }
    free(p->comment);
    p->comment = comment;
    
    /* read seed timecourse either from stdin or a supplied file */
    FILE *regressorStream = stdin;
    
    if ( p->regressorpath ) {
        regressorStream = fopen(p->regressorpath, "r");
    }
    
    p->regressor = rsReadRegressorFromStream(regressorStream, &p->nRegressorValues);

    if ( p->regressorpath ) {
        fclose(regressorStream);
    }
    
    /* open input file */
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }
    
    if ( p->verbose ) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "Regressor file: %s\n", (p->regressorpath==NULL ? "stdin" : p->regressorpath));
        fprintf(stdout, "Mask file: %s\n", p->maskpath);
        fprintf(stdout, "Seed length: %u\n", p->nRegressorValues);
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
        if (p->delay != 0) {
            fprintf(stdout, "Seed delay: %u\n", p->delay);
        }
        if (p->conversionMode == RSTOOLS_CORRELATION_CONVERSION_NONE) {
            fprintf(stdout, "Conversion: correlation coefficient\n");
        } else if (p->conversionMode == RSTOOLS_CORRELATION_CONVERSION_T) {
            fprintf(stdout, "Conversion: T-values\n");
        } else {
            fprintf(stdout, "Conversion: z-score\n");
        }
    }

    /* create correlation file (output) */
    p->output = rsCloneNiftiFile(p->outputpath, p->input, RSNIFTI_OPEN_ALLOC, 1);
    
    if ( ! p->output->readable ) {
        fprintf(stderr, "\nError: The nifti file containing the correlation output (%s) could not be created.\n", p->outputpath);
        return;
    }
    
    if (p->conversionMode == RSTOOLS_CORRELATION_CONVERSION_NONE) {
        FslSetIntent(p->output->fslio, NIFTI_INTENT_CORREL, p->input->vDim-2, 0, 0);
    } else if (p->conversionMode == RSTOOLS_CORRELATION_CONVERSION_T) {
        FslSetIntent(p->output->fslio, NIFTI_INTENT_TTEST, p->input->vDim-2, 0, 0);
    } else {
        FslSetIntent(p->output->fslio, NIFTI_INTENT_ZSCORE, 0, 0, 0);
    }
    
    /* prepare buffer */
    p->correlation = d3matrix(p->input->zDim-1, p->input->yDim-1, p->input->xDim-1);
    
    /* load mask */
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

void rsCorrelationRun(rsCorrelationParameters *p)
{   
    p->parametersValid = FALSE;
    
    short x,y,z,processedSlices = 0;
    double *timecourse;
    
    omp_lock_t updateProgressLock;
    omp_init_lock(&updateProgressLock);

    /* Iterate over all voxels for which the correlation is to be computed */
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(y,x,timecourse) shared(processedSlices)
    {
        #pragma omp for schedule(guided)
        for (short z=0; z<p->input->zDim; z++) {
            for (short y=0; y<p->input->yDim; y++) {
                for (short x=0; x<p->input->xDim; x++) {
                    
                    /* If it's not in the mask skip it to improve the performance */
                    if (p->mask != NULL && p->mask[z][y][x] < 0.1) {
                        
                        /* set the value in the correlation file to NaN so that it is skipped in later processing steps */
                        p->correlation[z][y][x] = log(-1.0);
                        
                        continue;
                    }
                    
                    /* read out timecourse */
                    double *fullTimecourse = malloc(p->input->vDim*sizeof(double));
                    rsExtractTimecourseFromRSNiftiFileBuffer(p->input, fullTimecourse, rsMakePoint3D(x, y, z));
                    
                    /* add the defined delay to the regressor and adjust the timecourse */
                    double *regressor = &p->regressor[(short)fmax(0, -1*p->delay)];
                    timecourse = &fullTimecourse[(short)fmax(0, p->delay)];
                    const size_t regressorLength = p->nRegressorValues - fabs(p->delay);
                    
                    /* compute correlation */
                    if ( p->monteCarloRepetitions > 0 ) {
                        p->correlation[z][y][x] = rsMonteCarloZCorrelation(timecourse, regressor, regressorLength, p->monteCarloRepetitions, p->monteCarloSampleSize);
                    } else if ( p->conversionMode == RSTOOLS_CORRELATION_CONVERSION_Z ) {
                        p->correlation[z][y][x] = rsZCorrelation(timecourse, regressor, regressorLength);
                    } else if ( p->conversionMode == RSTOOLS_CORRELATION_CONVERSION_NONE ) {
                        p->correlation[z][y][x] = rsCorrelation(timecourse, regressor, regressorLength);
                    } else if ( p->conversionMode == RSTOOLS_CORRELATION_CONVERSION_T ) {
                        p->correlation[z][y][x] = rsTCorrelation(timecourse, regressor, regressorLength);
                    }
                    
                    free(fullTimecourse);
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
    
    /* Write correlation file */
    rsCorrelationWriteCorrelationFile(p);

    p->parametersValid = TRUE;
}

void rsCorrelationWriteCorrelationFile(rsCorrelationParameters* p)
{
    
    /* Write correlation file */
    if (p->verbose) fprintf(stdout, "Writing correlation file\n");
    rsWriteNiftiHeader(p->output->fslio, p->comment);

    size_t buffsize = rsGetBufferSize(p->output->xDim, p->output->yDim, p->output->zDim, 1, p->output->dt);
    void *correlationBuffer;
    correlationBuffer = rsMalloc(buffsize);
    
    convertScaledDoubleToBuffer(
        p->output->dt,
        correlationBuffer,
        &p->correlation[0][0][0],
        p->output->slope,
        p->output->inter,
        p->output->xDim,
        p->output->yDim,
        p->output->zDim
    );
    
    FslWriteVolumes(p->output->fslio, correlationBuffer, 1);
    
    free(correlationBuffer);
}

void rsCorrelationDestroy(rsCorrelationParameters* p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFileAndFree(p->input);
        p->input = NULL;
    }
    
    if ( p->output != NULL ) {
        p->output->data = NULL;
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }
    
    if ( p->maskpath != NULL ) {
        free(p->mask);
        p->mask = NULL;
    }

    rsCorrelationFreeParams(p);
}

double *rsReadRegressorFromStream(FILE *stream, unsigned int *nValues)
{
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    double value;
    unsigned int nBuffer = 10;
    *nValues = 0;
    double *regressor = (double*)rsMalloc(nBuffer * sizeof(double));
    
    while ((read = getline(&line, &len, stream)) != -1) {
        value = atof(line);
        regressor[*nValues] = value;
        *nValues = *nValues + 1;
        
        // Check if we're running out of memory and extend the array if necessary
        if ( *nValues + 1 >= nBuffer ) {
            nBuffer = nBuffer + 10;
            double* tmpRegressor = realloc(regressor, nBuffer * sizeof(double));
            if (tmpRegressor) {
                regressor = tmpRegressor;
            } else {
                fprintf(stderr, "Could not allocate enough memory to save the regressor.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    
    if (line) free(line);
    
    return regressor;
}
