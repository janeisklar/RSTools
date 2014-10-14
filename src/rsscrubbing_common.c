#include "rsscrubbing_common.h"
#include "maths/utils.h"
#include "utils/rsio.h"

void rsScrubbingInit(rsScrubbingParameters *p)
{
    p->parametersValid = FALSE;

    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        (const char*)p->flaggedpath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->outputpath,
        RSIO_LASTFILE
    });
    
    if ( ! inputsReadable || ! outputsWritable ) {
        return;
    }
        
    /* open input file */
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }
    
    /* output the most important parameters to the user */
    if ( p->verbose ) {
        fprintf(stdout, "Input file:  %s\n", p->inputpath);
        fprintf(stdout, "Output file: %s\n", p->outputpath);
        fprintf(stdout, "Flagged frames file: %s\n", p->flaggedpath);
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
    }
    
    /* load realignment parameters */
    if ( ! rsLoadIndexVector(p->flaggedpath, &p->flaggedFrames, p->input->vDim) ) {
        fprintf(stderr, "\nError: The txt-file containing the flagged volume indices ('%s') could not be read!\n", p->flaggedpath);
        return;
    }

    p->parametersValid = TRUE;
}

void rsScrubbingRun(rsScrubbingParameters *p)
{
    p->parametersValid = FALSE;
    
    // count remaining frames
    int remainingFrames = 0;
    int framemap[p->input->vDim];

    for ( int t=0; t<p->input->vDim; t=t+1 ) {
        if ( p->flaggedFrames[t] == FALSE ) {
            framemap[remainingFrames] = t;
            remainingFrames = remainingFrames + 1;
        }
    }
    
    // stop here if the scrubbed file should not be saved
    if ( p->outputpath == NULL ) {
        p->parametersValid = TRUE;
        return;
    }

    // create output file
    p->output = rsCloneNiftiFile(p->outputpath, p->input, RSNIFTI_OPEN_ALLOC, remainingFrames);
    
    if ( ! p->output->readable ) {
        fprintf(stderr, "\nError: The nifti file containing the scrubbed output (%s) could not be created.\n", p->outputpath);
        return;
    }
    
    rsWriteNiftiHeader(p->output->fslio, p->callString);

    // prepare the output file's content
    int t;
    for ( t=0; t<remainingFrames; t++ ) {
        // extract a single volume for timepoint t from the buffer
        double ***data = d3matrix(p->input->zDim-1, p->input->yDim-1, p->input->xDim-1);
        rsExtractVolumeFromRSNiftiFileBuffer(p->input, data[0][0], framemap[t]);

        // write back to the buffer
        rsWriteVolumeToBuffer(p->output->dt, data[0][0], p->output->data, p->output->slope, p->output->inter, t, p->output->xDim, p->output->yDim, p->output->zDim);

        // free up memory
        free(data[0][0]); free(data[0]); free(data);
    }

    // write scrubbed file
    FslWriteVolumes(p->output->fslio, p->output->data, remainingFrames);

    p->parametersValid = TRUE;
}

void rsScrubbingDestroy(rsScrubbingParameters *p)
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

    rsScrubbingFreeParams(p);
}

/*
 * Load a vector of indices from a txt-file to a
 * BOOL vector containing TRUE if the index is contained
 * in the file and FALSE if it isn't
 */
BOOL rsLoadIndexVector(char *path, BOOL **indices, long length)
{
    BOOL *vector = (BOOL*)rsMalloc(sizeof(BOOL)*length);
    *indices = &vector[0];
    
    for (long i=0; i<length; i++) {
        vector[i] = FALSE;
    }
    
    unsigned int nValues;
    FILE *file = fopen(path, "r");
    double *content;
    content = rsReadRegressorFromStream(file, &nValues);
    fclose(file);
        
    for (unsigned int i=0; i<nValues; i++ ) {
        long value = (long)round(content[i]);

        if (value >= length) {
            return FALSE;
        }
        
        vector[value] = TRUE;
    }
    
    return TRUE;
}
