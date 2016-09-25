#include <nifti/headerinfo.h>
#include <nifti/rsniftiutils.h>
#include <externals/fslio/fslio.h>
#include <rstools/maths/geom.h>
#include "rsmaskborderdistance_common.h"
#include "rsmaskborderdistance_ui.h"
#include "utils/rsio.h"

void rsMaskBorderDistanceInit(rsMaskBorderDistanceParameters *p)
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

    if (p->verbose) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "Output file: %s\n", p->outputpath);
    }

    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);

    if (!p->input->readable) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }

    if (p->verbose) {
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
    }

    // Create output volume
    p->output = rsCloneNiftiFile(p->outputpath, p->input, RSNIFTI_OPEN_NONE, 1);
    
    if (!p->output->readable) {
        fprintf(stderr, "\nError: The nifti file containing the output (%s) could not be created.\n", p->outputpath);
        return;
    }

    // Change datatype
    FslSetDataType(p->output->fslio, DT_FLOAT32);
    p->output->pixtype = FslGetDataType(p->output->fslio, &(p->output->dt));

    // Allocate memory for the output volume
    p->output->data = rsMalloc(rsGetBufferSize(p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim, p->output->dt));

    if (p->output->data == NULL) {
        fprintf(stderr, "\nError: The memory for the output file could not be created.\n");
        return;
    }

    p->parametersValid = TRUE;
    return;
}

void rsMaskBorderDistanceRun(rsMaskBorderDistanceParameters *p)
{
    p->parametersValid = FALSE;

    // Prepare world matrix for coordinate conversion
    mat44 inputWorldMatrix = p->output->fslio->niftiptr->sform_code == NIFTI_XFORM_UNKNOWN
                             ? p->output->fslio->niftiptr->qto_xyz
                             : p->output->fslio->niftiptr->sto_xyz;

    // Prepare empty timecourse
    double emptybuffer[1];
    emptybuffer[0] = log(-1.0);

    // Read out mask
    double ***mask = d3matrix(p->input->zDim-1, p->input->yDim-1, p->input->xDim-1);
    rsExtractVolumeFromRSNiftiFileBuffer(p->input, mask[0][0], 0);

    if (p->verbose) {
        fprintf(stdout, "Go through voxels and determine whether they are in the mask or not\n");
    }

    // Count number of on and off-mask voxels and initialize those values with NaN
    short x, y, z;
    unsigned long nOffMaskVoxels = 0L;
    unsigned long nOnMaskVoxels = 0L;
    for (z=0; z<p->input->zDim; z++) {
        for (y = 0; y < p->input->yDim; y++) {
            for (x = 0; x < p->input->xDim; x++) {
                if (mask[z][y][x] < 0.1) {
                    nOffMaskVoxels += 1L;
                } else {
                    nOnMaskVoxels += 1L;
                }
            }
        }
    }

    // Pre-allocate list of off-mask voxel coordinates
    FloatPoint3D **offMaskPoints = (FloatPoint3D **)rsMalloc(sizeof(FloatPoint3D*) * nOffMaskVoxels);
    FloatPoint3D **onMaskPoints = (FloatPoint3D **)rsMalloc(sizeof(FloatPoint3D*) * nOnMaskVoxels);
    Point3D **onMaskPointsVX = (Point3D **)rsMalloc(sizeof(Point3D*) * nOnMaskVoxels);
    unsigned long offMaskVoxelIndex = 0L;
    unsigned long onMaskVoxelIndex = 0L;
    for (z=0; z<p->input->zDim; z++) {
        for (y = 0; y < p->input->yDim; y++) {
            for (x = 0; x < p->input->xDim; x++) {
                Point3D *point = rsMakePoint3D(x, y, z);
                FloatPoint3D *pointMM = rsMakeFloatPoint3D(0, 0, 0);

                // determine coordinate of this voxel
                FslGetMMCoord(inputWorldMatrix, x, y, z, &(pointMM->x), &(pointMM->y), &(pointMM->z));

                if (mask[z][y][x] < 0.1) {
                    offMaskPoints[offMaskVoxelIndex] = pointMM;
                    offMaskVoxelIndex += 1L;

                    // set the value in the filtered data to NaN
                    rsWriteTimecourseToRSNiftiFileBuffer(p->output, emptybuffer, point);

                    rsFree(point);
                } else {
                    onMaskPoints[onMaskVoxelIndex] = pointMM;
                    onMaskPointsVX[onMaskVoxelIndex] = point;
                    onMaskVoxelIndex += 1L;
                }
            }
        }
    }

    rsFree(mask[0][0]);
    rsFree(mask[0]);
    rsFree(mask);

    if (p->verbose) {
        fprintf(stdout, "Compute distances\n");
    }

    // Compute distance for every masked voxel
    double minDistance, currentDistance;

    #pragma omp parallel num_threads(rsGetThreadsNum()) private(minDistance,currentDistance,offMaskVoxelIndex,onMaskVoxelIndex) shared(nOffMaskVoxels,nOnMaskVoxels,onMaskPointsVX,onMaskPoints,offMaskPoints)
    {
        #pragma omp for schedule(guided)
        for (onMaskVoxelIndex = 0L; onMaskVoxelIndex < nOnMaskVoxels; onMaskVoxelIndex += 1L) {
            // find non-masked voxel with the highest distance
            minDistance = 1e10;

            for (offMaskVoxelIndex = 0L; offMaskVoxelIndex < nOffMaskVoxels; offMaskVoxelIndex += 1L) {
                // euclidean distance seems to be too slow :-(
                /*currentDistance = powf(
                    (
                        powf(onMaskPoints[onMaskVoxelIndex]->x - offMaskPoints[offMaskVoxelIndex]->x, 2.0) +
                        powf(onMaskPoints[onMaskVoxelIndex]->y - offMaskPoints[offMaskVoxelIndex]->y, 2.0) +
                        powf(onMaskPoints[onMaskVoxelIndex]->z - offMaskPoints[offMaskVoxelIndex]->z, 2.0)
                    ),
                    1.0/2.0
                );*/

                // so we use the manhattan distance instead
                currentDistance = fabsf(onMaskPoints[onMaskVoxelIndex]->x - offMaskPoints[offMaskVoxelIndex]->x) +
                                  fabsf(onMaskPoints[onMaskVoxelIndex]->y - offMaskPoints[offMaskVoxelIndex]->y) +
                                  fabsf(onMaskPoints[onMaskVoxelIndex]->z - offMaskPoints[offMaskVoxelIndex]->z);

                if (currentDistance < minDistance) {
                    minDistance = currentDistance;
                }
            }

            // write out filtered data to buffer
            rsWriteTimecourseToRSNiftiFileBuffer(p->output, &minDistance, onMaskPointsVX[onMaskVoxelIndex]);
        }
    }

    if (p->verbose) {
        fprintf(stdout, "Cleanup..\n");
    }

    for (offMaskVoxelIndex = 0L; offMaskVoxelIndex < nOffMaskVoxels; offMaskVoxelIndex++) {
        rsFree(offMaskPoints[offMaskVoxelIndex]);
    }
    rsFree(offMaskPoints);

    for (onMaskVoxelIndex = 0L; onMaskVoxelIndex < nOnMaskVoxels; onMaskVoxelIndex++) {
        rsFree(onMaskPoints[onMaskVoxelIndex]);
        rsFree(onMaskPointsVX[onMaskVoxelIndex]);
    }
    rsFree(onMaskPoints);
    rsFree(onMaskPointsVX);

    if ( p->verbose ) {
        fprintf(stdout, "Write out result to: %s\n", p->outputpath);
    }
    
    rsWriteNiftiHeader(p->output->fslio, p->callString);
    FslWriteVolumes(p->output->fslio, p->output->data, p->output->vDim);

    p->parametersValid = TRUE;
}

void rsMaskBorderDistanceDestroy(rsMaskBorderDistanceParameters *p)
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
    
    rsMaskBorderDistanceFreeParams(p);
}
