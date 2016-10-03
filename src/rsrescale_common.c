#include "rsrescale_common.h"
#include "maths/utils.h"
#include "utils/rsio.h"
#include "rsrescale_ui.h"
#include <math.h>

void rsResampleApplyInterpolation(double ***data_out, const rsRescaleParameters *p, double ***data_in, const short t, const short z);

void rsRescaleInit(rsRescaleParameters *p)
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
        
    /* open input file */
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }
    
    rsSetThreadsNum(p->threads);

    // adjust voxel dimensions
    float inputDim[] = {p->input->xDim, p->input->yDim, p->input->zDim};
    // determine dim sizes of the output
    p->newDim[0] = (short)ceil(((double)p->input->xDim) * p->scale[0]); // new x width
    p->newDim[1] = (short)ceil(((double)p->input->yDim) * p->scale[1]); // new y width
    p->newDim[2] = (short)ceil(((double)p->input->zDim) * p->scale[2]); // new z width

    // adjust voxel spacing
    FslGetVoxDim(p->input->fslio, &p->oldSpacing[0], &p->oldSpacing[1], &p->oldSpacing[2], &p->TR);
    for (short i=0; i<3; i++) {
        p->appliedScale[i] = ((double)p->newDim[i]) / ((double)inputDim[i]);
        p->newSpacing[i] = p->oldSpacing[i] / p->appliedScale[i];
    }

    // get coordinate matrix
    nifti_image *inputNifti = p->input->fslio->niftiptr;
    mat44 inputWorldMatrix = inputNifti->sform_code == NIFTI_XFORM_UNKNOWN
                             ? inputNifti->qto_xyz
                             : inputNifti->sto_xyz;

    // pre-compute inverse matrix for a faster computation of coordinates
    p->invInputWorldMatrix = nifti_mat44_inverse(inputWorldMatrix);

    // adjust coordinate matrices
    mat44 transformOriginMatrix;
    transformOriginMatrix.m[0][0] = 1;
    transformOriginMatrix.m[0][1] = 0;
    transformOriginMatrix.m[0][2] = 0;
    transformOriginMatrix.m[0][3] = inputWorldMatrix.m[0][3] * -1.0f;
    transformOriginMatrix.m[1][0] = 0;
    transformOriginMatrix.m[1][1] = 1;
    transformOriginMatrix.m[1][2] = 0;
    transformOriginMatrix.m[1][3] = inputWorldMatrix.m[1][3] * -1.0f;
    transformOriginMatrix.m[2][0] = 0;
    transformOriginMatrix.m[2][1] = 0;
    transformOriginMatrix.m[2][2] = 1;
    transformOriginMatrix.m[2][3] = inputWorldMatrix.m[2][3] * -1.0f;
    transformOriginMatrix.m[3][0] = 0;
    transformOriginMatrix.m[3][1] = 0;
    transformOriginMatrix.m[3][2] = 0;
    transformOriginMatrix.m[3][3] = 1;

    mat44 invTransformOriginMatrix;
    invTransformOriginMatrix = nifti_mat44_inverse(transformOriginMatrix);

    mat44 scalingMatrix;
    scalingMatrix.m[0][0] = 1.0 / p->appliedScale[0];
    scalingMatrix.m[0][1] = 0;
    scalingMatrix.m[0][2] = 0;
    scalingMatrix.m[0][3] = 0;
    scalingMatrix.m[1][0] = 0;
    scalingMatrix.m[1][1] = 1.0 / p->appliedScale[1];
    scalingMatrix.m[1][2] = 0;
    scalingMatrix.m[1][3] = 0;
    scalingMatrix.m[2][0] = 0;
    scalingMatrix.m[2][1] = 0;
    scalingMatrix.m[2][2] = 1.0 / p->appliedScale[0];
    scalingMatrix.m[2][3] = 0;
    scalingMatrix.m[3][0] = 0;
    scalingMatrix.m[3][1] = 0;
    scalingMatrix.m[3][2] = 0;
    scalingMatrix.m[3][3] = 1;

    // Transformations: translate to (0,0,0) origin, scale, then translate back to the original origin
    rsMat44MatrixMult(&scalingMatrix, &scalingMatrix, &transformOriginMatrix);
    rsMat44MatrixMult(&scalingMatrix, &invTransformOriginMatrix, &scalingMatrix);
    rsMat44MatrixMult(&p->outputWorldMatrix, &scalingMatrix, &inputWorldMatrix);

    // output the most important parameters to the user
    if ( p->verbose ) {
        fprintf(stdout, "Input file:  %s\n", p->inputpath);
        fprintf(stdout, "Output file: %s\n", p->outputpath);
		fprintf(stdout, "Input Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
		fprintf(stdout, "Output Dim: %d %d %d (%d Volumes)\n", p->newDim[0], p->newDim[1], p->newDim[2], p->input->vDim);
        fprintf(stdout, "Input voxel size: %.2f %.2f %.2f\n", p->oldSpacing[0], p->oldSpacing[1], p->oldSpacing[2]);
        fprintf(stdout, "Output voxel size: %.2f %.2f %.2f\n", p->newSpacing[0], p->newSpacing[1], p->newSpacing[2]);
        fprintf(stdout, "Applied scaling: %.4f %.4f %.4f (might slighty vary from the input as the scaling is applied such that it gives integer voxel dimensions)\n", p->appliedScale[0], p->appliedScale[1], p->appliedScale[2]);
        fprintf(stdout, "Old world matrix:\n");
        for (short i=0; i<4; i++) {
            fprintf(stdout, " %+.6f  %+.6f  %+.6f  %+.6f\n", inputWorldMatrix.m[i][0], inputWorldMatrix.m[i][1], inputWorldMatrix.m[i][2], inputWorldMatrix.m[i][3]);
        }
        fprintf(stdout, "\n");
        fprintf(stdout, "Scaling matrix:\n");
        for (short i=0; i<4; i++) {
            fprintf(stdout, " %+.6f  %+.6f  %+.6f  %+.6f\n", scalingMatrix.m[i][0], scalingMatrix.m[i][1], scalingMatrix.m[i][2], scalingMatrix.m[i][3]);
        }
        fprintf(stdout, "\n");
        fprintf(stdout, "New world matrix:\n");
        for (short i=0; i<4; i++) {
            fprintf(stdout, " %+.6f  %+.6f  %+.6f  %+.6f\n", p->outputWorldMatrix.m[i][0], p->outputWorldMatrix.m[i][1], p->outputWorldMatrix.m[i][2], p->outputWorldMatrix.m[i][3]);
        }
        fprintf(stdout, "\n");
    }
	
    p->parametersValid = TRUE;
}

void rsRescaleRun(rsRescaleParameters *p)
{
    p->parametersValid = FALSE;

	// create output file
	p->output = rsCloneNiftiFileWithNewDimensions(p->outputpath, p->input, RSNIFTI_OPEN_ALLOC, p->newDim[0], p->newDim[1], p->newDim[2], 0);

	if ( ! p->output->readable ) {
		fprintf(stderr, "\nError: The nifti file containing the output (%s) could not be created.\n", p->outputpath);
		return;
	}

    nifti_image *outputNifti = p->output->fslio->niftiptr;
    FslSetVoxDim(p->output->fslio, p->newSpacing[0], p->newSpacing[1], p->newSpacing[2], p->TR);
    FslSetRigidXform(p->output->fslio, outputNifti->qform_code, p->outputWorldMatrix);
    FslSetStdXform(p->output->fslio, outputNifti->sform_code, p->outputWorldMatrix);

	// write nifti header
	rsWriteNiftiHeader(p->output->fslio, p->callString);

    // prepare the output file's content
    int processedVolumes = 0;

    for (int t=0; t<p->input->vDim; t++ ) {

        // extract a single volume for timepoint t from the buffer
        double ***data_in = d3matrix(p->input->zDim - 1, p->input->yDim - 1, p->input->xDim - 1);
        rsExtractVolumeFromRSNiftiFileBuffer(p->input, data_in[0][0], t);

        double ***data_out = d3matrix(p->output->zDim-1, p->output->yDim-1, p->output->xDim-1);
        short z;

        #pragma omp parallel num_threads(rsGetThreadsNum()) private(z) shared(processedVolumes,data_in,data_out,t)
        {
            #pragma omp for schedule(guided, 1)
            for (z = 0; z < p->output->zDim; z++) {

                // interpolate volume at timepoint t
                rsResampleApplyInterpolation(data_out, p, data_in, t, z);
            }
	    }

        // write back to the buffer
        rsWriteVolumeToBuffer(p->output->dt, data_out[0][0], p->output->data, p->output->slope, p->output->inter, t, p->output->xDim, p->output->yDim, p->output->zDim);

        // free up memory
        free(data_in[0][0]);  free(data_in[0]);  free(data_in);
        free(data_out[0][0]); free(data_out[0]); free(data_out);

        // show progress
        if (p->progressCallback != NULL) {
            rsReportProgressEvent *event = (rsReportProgressEvent *) rsMalloc(sizeof(rsReportProgressEvent));
            event->run = processedVolumes;
            processedVolumes += 1;
            event->percentage = (double) processedVolumes * 100.0 / (double) p->input->vDim;
            rsReportProgressCallback_t cb = p->progressCallback->cb;
            void *data = p->progressCallback->data;
            cb(event, data);
            rsFree(event);
        } else if (p->verbose) {
            processedVolumes += 1;

            if (processedVolumes > 0 && processedVolumes % (short) (p->input->vDim >= 10 ? p->input->vDim / 10 : p->input->vDim) == 0) {
                fprintf(stdout, "..%.0f%%\n", ceil((float) processedVolumes * 100.0 / (float) p->input->vDim));
            }
        }
	}

    // write rescaled file
    FslWriteVolumes(p->output->fslio, p->output->data, p->input->vDim);

    p->parametersValid = TRUE;
}

void rsRescaleDestroy(rsRescaleParameters *p)
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

    rsRescaleFreeParams(p);
}

void rsResampleApplyInterpolation(double ***data_out, const rsRescaleParameters *p, double ***data_in, const short t, const short z)
{
    // convert pixdim to double (as required by rsInterpolation3DLanczosInterpolation()
    double pixdim[3] = {p->oldSpacing[0], p->oldSpacing[1], p->oldSpacing[2]};

    for (short y = 0; y < p->output->yDim; y++) {
        for (short x = 0; x < p->output->xDim; x++) {
            // get output coordinate in mm
            FloatPoint3D *outputPointInMM = rsMakeFloatPoint3D(0, 0, 0);
            FslGetMMCoord(p->outputWorldMatrix, x, y, z, &outputPointInMM->x, &outputPointInMM->y, &outputPointInMM->z);

            // get input coordinate in vx
            // Note: as we use the inverted world matrix the result will therefore
            // be in voxels rather than mm as the method name suggests
            FloatPoint3D *inputPointInVX = rsMakeFloatPoint3D(0, 0, 0);
            FslGetMMCoord(p->invInputWorldMatrix, outputPointInMM->x, outputPointInMM->y, outputPointInMM->z, &inputPointInVX->x, &inputPointInVX->y, &inputPointInVX->z);

            // interpolate current voxel
            if (p->linearInterpolation) {
                data_out[z][y][x] = rsTriLinearDistInterpolation(data_in, p->input->xDim, p->input->yDim, p->input->zDim, pixdim, inputPointInVX);
            } else {
                data_out[z][y][x] = rsInterpolation3DLanczosInterpolation(data_in, p->input->xDim, p->input->yDim, p->input->zDim, pixdim, inputPointInVX);
            }

            rsFree(outputPointInMM);
            rsFree(inputPointInVX);
        }
    }
}


