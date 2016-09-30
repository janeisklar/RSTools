#include <nifti/headerinfo.h>
#include <nifti/rsniftiutils.h>
#include <externals/fslio/fslio.h>
#include "maths/geom.h"
#include "rsmaskborderdistance_common.h"
#include "rssmoothing_common.h"
#include "rsmaskborderdistance_ui.h"
#include "utils/rsio.h"

void rsMarkNeighbouringVoxels(double ***result, double ***input, short xdim, short ydim, short zdim);

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

    // Get voxel spacing
    float xSpacing, ySpacing, zSpacing, tr;
    FslGetVoxDim(p->input->fslio, &xSpacing, &ySpacing, &zSpacing, &tr);


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

    // Binarize mask (>0)
    if (p->verbose) {
        fprintf(stdout, "Binarize mask\n");
    }

    short x, y, z;
    for (z=0; z<p->input->zDim; z++) {
        for (y = 0; y < p->input->yDim; y++) {
            for (x = 0; x < p->input->xDim; x++) {
                if (mask[z][y][x] < 0.1 || isnan(mask[z][y][x]) || isinf(mask[z][y][x])) {
                    mask[z][y][x] = 0.0;
                } else {
                    mask[z][y][x] = 1.0;
                }
            }
        }
    }

    // Create processing mask
    if (p->verbose) {
        fprintf(stdout, "Create processing mask by including mask neighbouring voxels only\n");
    }

    double ***processingMask = d3matrix(p->input->zDim-1, p->input->yDim-1, p->input->xDim-1);
    rsMarkNeighbouringVoxels(processingMask, mask, p->input->xDim, p->input->yDim, p->input->zDim);

    // Count number of on and off-mask voxels and initialize those values with NaN
    if (p->verbose) {
        fprintf(stdout, "Go through voxels and determine whether they are in the mask or not\n");
    }

    unsigned long nOffMaskVoxels = 0L;
    unsigned long nOnMaskVoxels = 0L;
    for (z=0; z<p->input->zDim; z++) {
        for (y = 0; y < p->input->yDim; y++) {
            for (x = 0; x < p->input->xDim; x++) {
                if (mask[z][y][x] > 0.1) {
                    nOnMaskVoxels += 1L;
                } else if (processingMask[z][y][x] > 0.3) {
                    nOffMaskVoxels += 1L;
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

                if (mask[z][y][x] > 0.1) {
                    onMaskPoints[onMaskVoxelIndex] = pointMM;
                    onMaskPointsVX[onMaskVoxelIndex] = point;
                    onMaskVoxelIndex += 1L;
                    continue;
                } else if (processingMask[z][y][x] > 0.3) {
                    offMaskPoints[offMaskVoxelIndex] = pointMM;
                    offMaskVoxelIndex += 1L;

                    // set the value in the filtered data to NaN
                    rsWriteTimecourseToRSNiftiFileBuffer(p->output, emptybuffer, point);
                }

                rsFree(point);
            }
        }
    }

    rsFree(mask[0][0]);           rsFree(mask[0]);           rsFree(mask);
    rsFree(processingMask[0][0]); rsFree(processingMask[0]); rsFree(processingMask);

    if (p->verbose) {
        fprintf(stdout, "Found %lu voxels in the mask and %lu outside of the mask\n", nOnMaskVoxels, nOffMaskVoxels);
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
                currentDistance = powf(
                    (
                        powf(onMaskPoints[onMaskVoxelIndex]->x - offMaskPoints[offMaskVoxelIndex]->x, 2.0) +
                        powf(onMaskPoints[onMaskVoxelIndex]->y - offMaskPoints[offMaskVoxelIndex]->y, 2.0) +
                        powf(onMaskPoints[onMaskVoxelIndex]->z - offMaskPoints[offMaskVoxelIndex]->z, 2.0)
                    ),
                    1.0/2.0
                );

                /*
                // so we use the manhattan distance instead
                currentDistance = fabsf(onMaskPoints[onMaskVoxelIndex]->x - offMaskPoints[offMaskVoxelIndex]->x) +
                                  fabsf(onMaskPoints[onMaskVoxelIndex]->y - offMaskPoints[offMaskVoxelIndex]->y) +
                                  fabsf(onMaskPoints[onMaskVoxelIndex]->z - offMaskPoints[offMaskVoxelIndex]->z);
                */

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

void rsMarkNeighbouringVoxels(double ***result, double ***input, short xdim, short ydim, short zdim) {
    for (short z=0; z<zdim; z++) {
        for (short y=0; y<ydim; y++) {
            for (short x=0; x<xdim; x++) {
                BOOL include = FALSE;
                for (short mz=-1; mz<2 && include == FALSE; mz++) {
                    const short z2 = z + mz;
                    for (short my=-1; my<2 && include == FALSE; my++) {
                        const short y2 = y + my;
                        for (short mx=-1; mx<2; mx++) {
                            const short x2 = x + mx;
                            if (x2>=0 && x2<xdim && y2>=0 && y2<ydim && z2>=0 && z2<zdim) {
                                if (input[z2][y2][x2] > 0.5) {
                                    include = TRUE;
                                    break;
                                }
                            }
                        }
                    }
                }
                result[z][y][x] = include ? 1.0 : 0.0;
            }
        }
    }
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
