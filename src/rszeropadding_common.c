#include <src/nifti/headerinfo.h>
#include <src/nifti/rsniftiutils.h>
#include <nifti1_io.h>
#include "nifti/rsniftiutils.h"
#include "rszeropadding_common.h"
#include "utils/rsio.h"
#include "rszeropadding_ui.h"

BOOL rsZeropaddingIsPaddedPoint(rsZeropaddingParameters *p, Point3D *point);
mat44 rsZeropaddingTranslateWorldMatrix(rsZeropaddingParameters *p, mat44 mat);

void rsZeropaddingInit(rsZeropaddingParameters *p)
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

    // compute size of the output dimensions
    p->newDim[0] = p->input->xDim + p->padding[0] + p->padding[1];
    p->newDim[1] = p->input->yDim + p->padding[2] + p->padding[3];
    p->newDim[2] = p->input->zDim + p->padding[4] + p->padding[5];

    if (p->newDim[0] <= 0) {
        fprintf(stderr, "\nError: The output nifti would have a negative size in x direction!\n");
        return;
    }

    if (p->newDim[1] <= 0) {
        fprintf(stderr, "\nError: The output nifti would have a negative size in y direction!\n");
        return;
    }

    if (p->newDim[2] <= 0) {
        fprintf(stderr, "\nError: The output nifti would have a negative size in z direction!\n");
        return;
    }
	
    /* output the most important parameters to the user */
    if ( p->verbose ) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "Output file: %s\n", p->outputpath);
        fprintf(stdout, "Input Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
        fprintf(stdout, "Output Dim: %d %d %d (%d Volumes)\n", p->newDim[0], p->newDim[1], p->newDim[2], p->input->vDim);
        if (p->mirroredPadding) {
            fprintf(stdout, "Padding will be performed by mirroring the edges\n");
        } else {
            fprintf(stdout, "Value used for padding: %.2f\n", p->paddingValue);
        }
    }
	
    p->parametersValid = TRUE;
}

void rsZeropaddingRun(rsZeropaddingParameters *p)
{
    p->parametersValid = FALSE;

	// get voxel spacing
	float xSpacing, ySpacing, zSpacing, tr;
	FslGetVoxDim(p->input->fslio, &xSpacing, &ySpacing, &zSpacing, &tr);

    // create output file
    p->output = rsCloneNiftiFileWithNewDimensions(p->outputpath, p->input, RSNIFTI_OPEN_ALLOC, p->newDim[0], p->newDim[1], p->newDim[2], 0);

    if ( ! p->output->readable ) {
        fprintf(stderr, "\nError: The nifti file containing the output (%s) could not be created.\n", p->outputpath);
        return;
    }

    // convert world matrix
    nifti_image *outputNifti = p->output->fslio->niftiptr;
    outputNifti->qto_xyz = rsZeropaddingTranslateWorldMatrix(p, outputNifti->qto_xyz);
    outputNifti->sto_xyz = rsZeropaddingTranslateWorldMatrix(p, outputNifti->sto_xyz);

    FslSetRigidXform(p->output->fslio, outputNifti->qform_code, outputNifti->qto_xyz);
    FslSetStdXform(p->output->fslio, outputNifti->sform_code, outputNifti->sto_xyz);

    // write nifti header
	rsWriteNiftiHeader(p->output->fslio, p->callString);

    // prepare padding content buffer
    const size_t paddingBufferSize = rsGetBufferSize(1, 1, 1, p->output->vDim, p->output->dt);
    void *paddingBuffer = rsMalloc(paddingBufferSize);
    double *paddingValues = rsMalloc(p->output->vDim * sizeof(double));
    for (short i=0; i<p->output->vDim; i++) {
        paddingValues[i] = p->paddingValue;
    }
    convertScaledDoubleToBuffer(p->input->dt, paddingBuffer, paddingValues, p->output->slope, p->output->inter, 1, 1, 1);

    // prepare the output file's content
    unsigned short x,y,z;
    Point3D *pointIn;
    Point3D *pointOut;
    SignedPoint3D *signedPointIn;

    for (z=0; z<p->output->zDim; z++) {
        for (y = 0; y < p->output->yDim; y++) {
            for (x = 0; x < p->output->xDim; x++) {

                pointOut = rsMakePoint3D(x, y, z);

                if (rsZeropaddingIsPaddedPoint(p, pointOut)) {
                    if (p->mirroredPadding) {
                        // fill with mirrored version of the input
                        signedPointIn = rsMakeSignedPoint3D(x-p->padding[0], y-p->padding[2], z-p->padding[4]);
                        if (signedPointIn->x < 0 || signedPointIn->x >= (p->input->xDim - 1)) {
                            signedPointIn->x = signedPointIn->x < 0 ? -1 - signedPointIn->x : 2*p->input->xDim - signedPointIn->x - 1;
                        }
                        if (signedPointIn->y < 0 || signedPointIn->y >= (p->input->yDim - 1)) {
                            signedPointIn->y = signedPointIn->y < 0 ? -1 - signedPointIn->y : 2*p->input->yDim - signedPointIn->y - 1;
                        }
                        if (signedPointIn->z < 0 || signedPointIn->z >= (p->input->zDim - 1)) {
                            signedPointIn->z = signedPointIn->z < 0 ? -1 - signedPointIn->z : 2*p->input->zDim - signedPointIn->z - 1;
                        }
                        pointIn = rsMakePoint3D(signedPointIn->x, signedPointIn->y, signedPointIn->z);
                        rsCopyTimecourseFromInBufferToOutBuffer(
                            p->input->dt,     // datatype
                            p->output->data,  // output specification
                            pointOut,
                            p->output->xDim,
                            p->output->yDim,
                            p->output->zDim,
                            p->output->vDim,
                            p->input->data,   // input specification
                            pointIn,
                            p->input->xDim,
                            p->input->yDim,
                            p->input->zDim
                        );
                        rsFree(signedPointIn);
                    } else {
                        // fill with padding values
                        pointIn = rsMakePoint3D(0, 0, 0);
                        rsCopyTimecourseFromInBufferToOutBuffer(
                            p->input->dt,     // datatype
                            p->output->data,  // output specification
                            pointOut,
                            p->output->xDim,
                            p->output->yDim,
                            p->output->zDim,
                            p->output->vDim,
                            paddingBuffer,    // input specification (the padding buffer)
                            pointIn,
                            1,
                            1,
                            1
                        );
                    }
                } else {
                    // fill with timecourse from the corresponding point in the input volume
                    pointIn = rsMakePoint3D(x-p->padding[0], y-p->padding[2], z-p->padding[4]);
                    rsCopyTimecourseFromInBufferToOutBuffer(
                        p->input->dt,     // datatype
                        p->output->data,  // output specification
                        pointOut,
                        p->output->xDim,
                        p->output->yDim,
                        p->output->zDim,
                        p->output->vDim,
                        p->input->data,   // input specification
                        pointIn,
                        p->input->xDim,
                        p->input->yDim,
                        p->input->zDim
                    );
                }

                rsFree(pointIn);
                rsFree(pointOut);
            }
        }
    }

    // write transformed file
    FslWriteVolumes(p->output->fslio, p->output->data, p->output->vDim);

    p->parametersValid = TRUE;
}

void rsZeropaddingDestroy(rsZeropaddingParameters *p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFileAndFree(p->input);
        p->input = NULL;
    }
    
    if ( p->output != NULL ) {
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }

    rsZeropaddingFreeParams(p);
}

/**
 * Check if the point still lies outside of the input volume
 * NOTE: will not check for a negative upper padding as this function won't
 *       be called in that case anyway. padding[1], padding[2] and padding[3]
 *       are therefore ignored.
 */
BOOL rsZeropaddingIsPaddedPoint(rsZeropaddingParameters *p, Point3D *point)
{
    if ((p->padding[0] > 0 && point->x < p->padding[0]) || (point->x >= (p->input->xDim + p->padding[0]))) {
        return TRUE;
    }

    if ((p->padding[2] > 0 && point->y < p->padding[2]) || (point->y >= (p->input->yDim + p->padding[2]))) {
        return TRUE;
    }

    if ((p->padding[4] > 0 && point->z < p->padding[4]) || (point->z >= (p->input->zDim + p->padding[4]))) {
        return TRUE;
    }

    return FALSE;
}

mat44 rsZeropaddingTranslateWorldMatrix(rsZeropaddingParameters *p, mat44 mat)
{
    mat44 result = mat;

    const float dx = p->padding[0];
    const float dy = p->padding[0];
    const float dz = p->padding[0];

    // fix translational params
    result.m[0][3] -= mat.m[0][0] * dx + mat.m[0][1] * dy + mat.m[0][2] * dz;
    result.m[1][3] -= mat.m[1][0] * dx + mat.m[1][1] * dy + mat.m[1][2] * dz;
    result.m[2][3] -= mat.m[2][0] * dx + mat.m[2][1] * dy + mat.m[2][2] * dz;
    result.m[3][3] -= mat.m[3][0] * dx + mat.m[3][1] * dy + mat.m[3][2] * dz;

    return result;
}