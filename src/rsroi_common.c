#include "rsroi_common.h"

void rsRoiInit(rsRoiParameters *p)
{
    p->parametersValid = FALSE;

    // check if the required arguments have been provided
    if ( p->inputpath == NULL ) {
        fprintf(stderr, "No input volume specified!\n");
        return;
    }

    if ( p->maskpath == NULL ) {
        fprintf(stderr, "A  path for the resulting mask must be specified!\n");
        return;
    }

    if ( p->center->x < -9999.0 && p->nSamples < 0 ) {
        fprintf(stderr, "ROI center needs to be specified!(-center)!\n");
        return;
    }

    if ( p->sphereradius <= 0 && p->cubeDim->x < 0 && p->nSamples < 0 ) {
        fprintf(stderr, "ROI sphere radius or cube dimensions needs to be specified!(--sphere, --cube)!\n");
        return;
    }

    p->inputEqualsOutput = ! strcmp(p->inputpath, p->maskpath);

    rsSetThreadsNum(1);

    // open input file
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_NONE);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }

    p->input->data = rsMalloc(rsGetBufferSize(p->input->xDim, p->input->yDim, p->input->zDim, 1, p->input->dt));
    p->input->vDim = 1;
    FslReadVolumes(p->input->fslio, p->input->data, p->input->vDim);
    rsCloseNiftiFile(p->input, TRUE);

    // prepare mask file
    p->mask = rsCloneNiftiFile(p->maskpath, p->input, RSNIFTI_CLONE_MEMORY, 1);

    // output the most important parameters to the user
    if ( p->verbose ) {
        char *unit = "mm";
        if ( p->useImageSpaceCoordinates ) {
            unit = "vx";
        }
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "Mask file: %s\n", p->maskpath);
        if ( p->sphereradius > 0 || p->cubeDim->x > -9999.9 )
        fprintf(stdout, "Center: %.2f%s %.2f%s %.2f%s\n", p->center->x, unit, p->center->y, unit, p->center->z, unit);
        if ( p->sphereradius > 0) {
            fprintf(stdout, "Sphere radius: %.4f%s\n", p->sphereradius, unit);
        }
        if ( p->cubeDim->x > -9999.9 ) {
            fprintf(stdout, "Cube: %.2f%s %.2f%s %.2f%s\n", p->cubeDim->x, unit, p->cubeDim->y, unit, p->cubeDim->z, unit);
        }
    }

    p->parametersValid = TRUE;
}

void rsRoiRun(rsRoiParameters *p)
{
    /* Get sForm/qForm for MM<->Voxel space conversion */
    int sform_code, qform_code;
    mat44 sform44, qform44, stdmat44;

    sform_code = FslGetStdXform(p->mask->fslio, &sform44);
    qform_code = FslGetRigidXform(p->mask->fslio, &qform44);

    if (sform_code!=NIFTI_XFORM_UNKNOWN) {
        stdmat44 = sform44;
    } else if (qform_code!=NIFTI_XFORM_UNKNOWN) {
        stdmat44 = qform44;
    }

    // Iterate over all voxels
    for (short z=0; z<p->input->zDim; z=z+1) {
        for (short y=0; y<p->input->yDim; y=y+1) {
            for (short x=0; x<p->input->xDim; x=x+1) {

                float mmx = x, mmy = y, mmz = z;

                if ( ! p->useImageSpaceCoordinates ) {
                    // Convert current coordinate to MM space
                    FslGetMMCoord(stdmat44, x, y, z, &mmx, &mmy, &mmz);
                }

                FloatPoint3D *point = rsMakeFloatPoint3D(mmx, mmy, mmz);

                double value;

                if ( p->sphereradius > 0 && rsVoxelInSphere(point, p->center, p->sphereradius) ) {
                    value = p->roiValue;
                } else if(p->cubeDim->x > -9999.9 && rsVoxelInCube(point, p->center, p->cubeDim)) {
                    value = p->roiValue;
                } else if ( ! p->keepVolume ) {
                    value = 0.0;
                } else {
                    rsFree(point);
                    continue;
                }

                Point3D *maskPoint = rsMakePoint3D(x, y, z);

                rsWriteTimecourseToRSNiftiFileBuffer(p->mask, &value, maskPoint);

                rsFree(point);
                rsFree(maskPoint);
            }
        }
    }

    // If sampling is requested sample nSamples voxels from the mask
    if ( p->nSamples > 0 ) {
        if ( p->verbose ) fprintf(stdout, "Sampling %ld voxels\n", p->nSamples);

        // Create a new empty (misuse the input volume for this matter)
        const double v = 0.0;
        for (short z=0; z<p->input->zDim; z=z+1) {
            for (short y=0; y<p->input->yDim; y=y+1) {
                for (short x=0; x<p->input->xDim; x=x+1) {
                    Point3D *maskPoint = rsMakePoint3D(x, y, z);
                    rsWriteTimecourseToRSNiftiFileBuffer(p->input, &v, maskPoint);
                    rsFree(maskPoint);
                }
            }
        }

        // Randomly take sample points and copy them from the original mask
        gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
        double oldValue, newValue;

        while (p->nSamples > 0) {
            unsigned short x = gsl_rng_uniform_int(r, p->input->xDim);
            unsigned short y = gsl_rng_uniform_int(r, p->input->yDim);
            unsigned short z = gsl_rng_uniform_int(r, p->input->zDim);

            Point3D *maskPoint = rsMakePoint3D(x, y, z);

            rsExtractTimecourseFromRSNiftiFileBuffer(p->mask,  &oldValue, maskPoint);
            rsExtractTimecourseFromRSNiftiFileBuffer(p->input, &newValue, maskPoint);

            if (oldValue <= 0.5) {
                rsFree(maskPoint);
                continue;
            }

            if (newValue > 0.5) {
                rsFree(maskPoint);
                continue;
            }

            rsWriteTimecourseToRSNiftiFileBuffer(p->input, &p->roiValue, maskPoint);
            p->nSamples = p->nSamples - 1L;

            rsFree(maskPoint);
        }

        gsl_rng_free(r);

        void *data = p->mask->data;
        p->mask->data = p->input->data;
        p->input->data = data;
    }

    // write file
    rsWriteNiftiHeader(p->mask->fslio, p->callString);
    FslWriteVolumes(p->mask->fslio, p->mask->data, p->mask->vDim);
}

void rsRoiDestroy(rsRoiParameters *p)
{
    if ( p->input != NULL ) {
        if ( p->input->data != NULL ) {
            rsFree(p->input->data);
        }
        rsFreeNiftiFile(p->input);
        p->input = NULL;
    }

    if ( p->mask != NULL ) {
        p->mask->data = NULL;
        rsCloseNiftiFileAndFree(p->mask);
        p->mask = NULL;
    }

    rsRoiFreeParams(p);
}
