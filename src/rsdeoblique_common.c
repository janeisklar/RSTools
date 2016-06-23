#include <src/nifti/headerinfo.h>
#include <src/nifti/rsniftiutils.h>
#include <nifti1_io.h>
#include <externals/fslio/fslio.h>
#include "nifti/rsniftiutils.h"
#include "rsdeoblique_common.h"
#include "utils/rsio.h"
#include "rsdeoblique_ui.h"
#include "math.h"

void rsDeobliqueWorldMatrix(mat44 *output, mat44 *transform, short newDims[3], const mat44 *input, const rsDeobliqueParameters *p);

void rsDeobliqueInit(rsDeobliqueParameters *p)
{
    p->parametersValid = FALSE;

    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->outputpath,
        (const char*)p->transformationpath,
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

    // determine voxel size of the output
    float pixDims[3];
    FslGetVoxDim(p->input->fslio, &pixDims[0], &pixDims[1], &pixDims[2], &p->TR);
    p->pixDimIn[0] = pixDims[0];
    p->pixDimIn[1] = pixDims[1];
    p->pixDimIn[2] = pixDims[2];
    if (p->retainVoxelSizes) {
        p->pixDimOut[0] = p->pixDimIn[0];
        p->pixDimOut[1] = p->pixDimIn[1];
        p->pixDimOut[2] = p->pixDimIn[2];
    } else {
        double pixDimOut = fmin(fmin(p->pixDimIn[0], p->pixDimIn[1]), p->pixDimIn[2]);
        p->pixDimOut[0] = pixDimOut;
        p->pixDimOut[1] = pixDimOut;
        p->pixDimOut[2] = pixDimOut;
    }

    if (p->pixDimOut[0] <= 0.001 || p->pixDimOut[1] <= 0.001 || p->pixDimOut[2] <= 0.001) {
        fprintf(stderr, "\nError: The voxel size of the input nifti must not be 0mm in any direction (x,y,z)!\n");
        return;
    }

    // determine dim sizes of the output
    p->dimsOut[0] = (short)ceil(((float)p->input->xDim) * (p->pixDimIn[0] / p->pixDimOut[0])); // new x width
    p->dimsOut[1] = (short)ceil(((float)p->input->yDim) * (p->pixDimIn[1] / p->pixDimOut[1])); // new y width
    p->dimsOut[2] = (short)ceil(((float)p->input->zDim) * (p->pixDimIn[2] / p->pixDimOut[2])); // new z width

    /* output the most important parameters to the user */
    if ( p->verbose ) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "Output file: %s\n", p->outputpath);
        fprintf(stdout, "Input voxel size: %.2f %.2f %.2f\n", p->pixDimIn[0], p->pixDimIn[1], p->pixDimIn[2]);
        fprintf(stdout, "Output voxel size: %.2f %.2f %.2f\n", p->pixDimOut[0], p->pixDimOut[1], p->pixDimOut[2]);
        fprintf(stdout, "Input Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
    }
	
    p->parametersValid = TRUE;
}

void rsDeobliqueRun(rsDeobliqueParameters *p)
{
    p->parametersValid = FALSE;

    // convert world matrix
    nifti_image *inputNifti = p->input->fslio->niftiptr;
    mat44 inputWorldMatrix = inputNifti->sform_code == NIFTI_XFORM_UNKNOWN
                             ? inputNifti->qto_xyz
                             : inputNifti->sto_xyz;
    mat44 outputWorldMatrix;
    mat44 transformMatrix;
    rsDeobliqueWorldMatrix(&outputWorldMatrix, &transformMatrix, &p->dimsOut[0], &inputWorldMatrix, p);

    if (p->verbose) {
        fprintf(stdout, "Output Dim: %d %d %d (%d Volumes)\n", p->dimsOut[0], p->dimsOut[1], p->dimsOut[2], p->input->vDim);
        fprintf(stdout, "Old world matrix:\n");
        for (short i=0; i<4; i++) {
            fprintf(stdout, " %+.6f  %+.6f  %+.6f  %+.6f\n", inputWorldMatrix.m[i][0], inputWorldMatrix.m[i][1], inputWorldMatrix.m[i][2], inputWorldMatrix.m[i][3]);
        }
        fprintf(stdout, "\n");
        fprintf(stdout, "Transform matrix:\n");
        for (short i=0; i<4; i++) {
            fprintf(stdout, " %+.6f  %+.6f  %+.6f  %+.6f\n", transformMatrix.m[i][0], transformMatrix.m[i][1], transformMatrix.m[i][2], transformMatrix.m[i][3]);
        }
        fprintf(stdout, "\n");
        fprintf(stdout, "New world matrix:\n");
        for (short i=0; i<4; i++) {
            fprintf(stdout, " %+.6f  %+.6f  %+.6f  %+.6f\n", outputWorldMatrix.m[i][0], outputWorldMatrix.m[i][1], outputWorldMatrix.m[i][2], outputWorldMatrix.m[i][3]);
        }
        fprintf(stdout, "\n");
    }

    // create output file
    p->output = rsCloneNiftiFileWithNewDimensions(p->outputpath, p->input, RSNIFTI_OPEN_ALLOC, p->dimsOut[0], p->dimsOut[1], p->dimsOut[2], 0);
    nifti_image *outputNifti = p->output->fslio->niftiptr;

    if ( ! p->output->readable ) {
        fprintf(stderr, "\nError: The nifti file containing the output (%s) could not be created.\n", p->outputpath);
        return;
    }

    FslSetVoxDim(p->output->fslio, p->pixDimOut[0], p->pixDimOut[1], p->pixDimOut[2], p->TR);
    FslSetRigidXform(p->output->fslio, outputNifti->qform_code, outputWorldMatrix);
    FslSetStdXform(p->output->fslio, outputNifti->sform_code, outputWorldMatrix);

    // pre-compute inverse matrix for a faster computation of coordinates
    mat44 invInputWorldMatrix;
    invInputWorldMatrix = nifti_mat44_inverse(inputWorldMatrix);

    // write nifti header
	rsWriteNiftiHeader(p->output->fslio, p->callString);

    // write out transform matrix
    if (p->transformationpath != NULL) {
        p->transform = fopen(p->transformationpath, "w+");

        if ( ! p->transform ) {
            fprintf(stderr, "\nError: The transformation file (%s) could not be created.\n", p->transformationpath);
            return;
        }

        fprintf(p->transform, "#Insight Transform File V1.0\n");
        fprintf(p->transform, "#Transform 0\n");
        fprintf(p->transform, "Transform: AffineTransform_double_3_3\n");
        fprintf(p->transform, "Parameters:");
        for (short i=0; i<3; i++) {
            fprintf(p->transform, " %.14f %.14f %.14f", transformMatrix.m[i][0], transformMatrix.m[i][1], transformMatrix.m[i][2]);
        }
        for (short i=0; i<3; i++) {
            fprintf(p->transform, " %.14f", transformMatrix.m[i][3]);
        }
        fprintf(p->transform, "\n");
        fprintf(p->transform, "FixedParameters: %.14f %.14f %.14f\n", 0.0, 0.0, 0.0);
    }

    int t;
    for ( t=0; t<p->input->vDim; t++ ) {

        // extract a single volume for timepoint t from the buffer
        double ***data_in = d3matrix(p->input->zDim-1, p->input->yDim-1, p->input->xDim-1);
        rsExtractVolumeFromRSNiftiFileBuffer(p->input, data_in[0][0], t);

        // interpolate in x direction
        double ***data_out = d3matrix(p->output->zDim-1, p->output->yDim-1, p->output->xDim-1);
        for (short z = 0; z < p->output->zDim; z++) {
            for (short y = 0; y < p->output->yDim; y++) {
                for (short x = 0; x < p->output->xDim; x++) {
                    // get output coordinate in mm
                    FloatPoint3D *outputPointInMM = rsMakeFloatPoint3D(0, 0, 0);
                    FslGetMMCoord(outputWorldMatrix, x, y, z, &outputPointInMM->x, &outputPointInMM->y, &outputPointInMM->z);

                    // get input coordinate in vx
                    // Note: as we use the inverted world matrix the result will therefore
                    // be in voxels rather than mm as the method name suggests
                    FloatPoint3D *inputPointInVX = rsMakeFloatPoint3D(0, 0, 0);
                    FslGetMMCoord(invInputWorldMatrix, outputPointInMM->x, outputPointInMM->y, outputPointInMM->z, &inputPointInVX->x, &inputPointInVX->y, &inputPointInVX->z);

                    // interpolate current voxel
                    //data_out[z][y][x] = rsTriLinearDistInterpolation(data_in, p->input->xDim, p->input->yDim, p->input->zDim, p->pixDimIn, inputPointInVX);
                    data_out[z][y][x] = rsInterpolation3DLanczosInterpolation(data_in, p->input->xDim, p->input->yDim, p->input->zDim, p->pixDimIn, inputPointInVX);
                    //data_out[z][y][x] = rsInterpolationTriLinearInterpolation(data_in, p->input->xDim, p->input->yDim, p->input->zDim, p->pixDimIn, inputPointInVX);

                    rsFree(outputPointInMM);
                    rsFree(inputPointInVX);
                }
            }
        }

        // write back to the buffer
        rsWriteVolumeToBuffer(p->output->dt, data_out[0][0], p->output->data, p->output->slope, p->output->inter, t, p->output->xDim, p->output->yDim, p->output->zDim);

        // free up memory
        free(data_in[0][0]);  free(data_in[0]);  free(data_in);
        free(data_out[0][0]); free(data_out[0]); free(data_out);
    }

    // write transformed file
    FslWriteVolumes(p->output->fslio, p->output->data, p->output->vDim);

    p->parametersValid = TRUE;
}

void rsDeobliqueDestroy(rsDeobliqueParameters *p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFileAndFree(p->input);
        p->input = NULL;
    }
    
    if ( p->output != NULL ) {
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }

    if ( p->transform != NULL ) {
        fclose(p->transform);
        p->transform = NULL;
    }

    rsDeobliqueFreeParams(p);
}

void rsDeobliqueWorldMatrix(mat44 *output, mat44 *transform, short newDims[3], const mat44 *input, const rsDeobliqueParameters *p)
{
    mat44 desiredOutput;

    // set desired output
    for (short i=0; i<4; i++) {
        for (short j=0; j<4; j++) {
            desiredOutput.m[i][j] = i==j ? p->pixDimOut[i] : 0;
        }
    }
    desiredOutput.m[0][3] = input->m[0][3];
    desiredOutput.m[1][3] = input->m[1][3];
    desiredOutput.m[2][3] = input->m[2][3];
    desiredOutput.m[3][3] = input->m[3][3];

    // compute inverse of input matrix
    mat44 invInput = nifti_mat44_inverse(*input);

    // compute required transformation to make it deoblique: transform = inv(input) * desiredOutput
    rsMat44MatrixMult(transform, &invInput, &desiredOutput);
    rsMat44MatrixMult(output, transform, input);

    // compute the boundary size necessary to fit the whole volume after rotation
    FloatPoint3D *inputBoundariesVX[8] = {
        rsMakeFloatPoint3D(0, 0, 0),
        rsMakeFloatPoint3D(p->input->xDim-1, 0, 0),
        rsMakeFloatPoint3D(0, p->input->yDim-1, 0),
        rsMakeFloatPoint3D(p->input->xDim-1, p->input->yDim-1, 0),
        rsMakeFloatPoint3D(0, 0, p->input->zDim-1),
        rsMakeFloatPoint3D(p->input->xDim-1, 0, p->input->zDim-1),
        rsMakeFloatPoint3D(0, p->input->yDim-1, p->input->zDim-1),
        rsMakeFloatPoint3D(p->input->xDim-1, p->input->yDim-1, p->input->zDim-1)
    };

    short min_x=0, max_x=0,
          min_y=0, max_y=0,
          min_z=0, max_z=0;

    for (short i=0; i<8; i++) {
        FloatPoint3D *inputBoundaryMM = rsMakeFloatPoint3D(0, 0, 0);
        FloatPoint3D *outputBoundaryVX = rsMakeFloatPoint3D(0, 0, 0);
        FslGetMMCoord(*input, inputBoundariesVX[i]->x, inputBoundariesVX[i]->y, inputBoundariesVX[i]->z, &inputBoundaryMM->x, &inputBoundaryMM->y, &inputBoundaryMM->z);
        FslGetVoxCoord(*output, inputBoundaryMM->x, inputBoundaryMM->y, inputBoundaryMM->z, &outputBoundaryVX->x, &outputBoundaryVX->y, &outputBoundaryVX->z);

        if (floor(outputBoundaryVX->x) < min_x)
            min_x = floor(outputBoundaryVX->x);
        if (ceil(outputBoundaryVX->x) > max_x)
            max_x = ceil(outputBoundaryVX->x);
        if (floor(outputBoundaryVX->y) < min_y)
            min_y = floor(outputBoundaryVX->y);
        if (ceil(outputBoundaryVX->y) > max_y)
            max_y = ceil(outputBoundaryVX->y);
        if (floor(outputBoundaryVX->z) < min_z)
            min_z = floor(outputBoundaryVX->z);
        if (ceil(outputBoundaryVX->z) > max_z)
            max_z = ceil(outputBoundaryVX->z);
    }

    short tmpDims[3] = {
        abs(max_x-min_x),
        abs(max_y-min_y),
        abs(max_z-min_z)
    };

    fprintf(stdout, "Padding output by %d %d %d voxels\n", tmpDims[0]-p->dimsOut[0], tmpDims[1]-p->dimsOut[1], tmpDims[2]-p->dimsOut[2]);
    newDims[0] = tmpDims[0];
    newDims[1] = tmpDims[1];
    newDims[2] = tmpDims[2];

    // compute the translational shift required due to padding
    float shiftInMM[3] = {0,0,0};

    FslGetMMCoord(*output, min_x, min_y, min_z, &shiftInMM[0], &shiftInMM[1], &shiftInMM[2]); // coordinate of shift
    shiftInMM[0] -= output->m[0][3]; // remove origin
    shiftInMM[1] -= output->m[1][3];
    shiftInMM[2] -= output->m[2][3];

    fprintf(stdout, "Shifting by %d %d %d voxels (%.2f %.2f %.2f mm) due to padding\n", abs(min_x), abs(min_y), abs(min_z), shiftInMM[0], shiftInMM[1], shiftInMM[2]);

    // create and apply translation transformation to account for the padding
    mat44 translationMatrix;
    for (short i=0; i<4; i++) {
        for (short j=0; j<4; j++) {
            translationMatrix.m[i][j] = i==j ? 1 : 0;
        }
    }
    translationMatrix.m[0][3] = shiftInMM[0];
    translationMatrix.m[1][3] = shiftInMM[1];
    translationMatrix.m[2][3] = shiftInMM[2];

    rsMat44MatrixMult(transform, &translationMatrix, transform);
    rsMat44MatrixMult(output, transform, input);
}