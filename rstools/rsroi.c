//
//  rsroi.c
//  rstools
//
//  Created by André Hoffmann on 6/7/13.
//
//
#include <stdio.h>
#include <nifti1.h>
#include <fslio.h>
#include "rsniftiutils.h"
#include "rsmathutils.h"

int show_help( void )
{
    printf(
       "rsroi: Given a 4D-Nifti that will be cloned this tool will create\n"
       "       a binary mask for the specified region(sphere or cube).\n"
       "\n"
    );
    
    printf(
       "basic usage:  rsroi -input <volume> -mask <volume> -sphere <radius> -center <x> <y> <z> [-verbose]\n"
       "\n"
    );
    
    printf(
       "options:\n"
    );
    
    printf(
       "   -help               : show this help\n"
    );
    
    printf(
       "   -input <volume>     : the volume from which the header, dimension and\n"
       "                         alignment infos will be taken. The resulting mask\n"
       "                         will thus be coregistered to this volume.\n"
    );
    
    printf(
       "   -sphere <radius>    : use this option if the ROI that is to be added is\n"
       "                         a sphere. Its radius should be given in mm.\n"
    );
    
    printf(
       "   -cube <x> <y> <z>   : use this option if the ROI that is to be added is\n"
       "                         a cube. The height/width/depth should be given in mm.\n"
    );
    
    printf(
       "   -center <x> <y> <z> : this option specifies the center of the ROI. \n"
    );
    
    printf(
       "   -mask <mask>        : a mask specifying the ROI for improved performance\n"
    );
    
    printf(
       "   -keepVolume         : keep values from the input volume. Use this option\n"
       "                         if you want to add ROIs to a mask that you specified\n"
       "                         as an input.\n"
    );
    
    printf(
        "  -randomsample <n>   : randomly sample <n> voxels from the file specified\n"
        "                        using -input(through -keepVolume) or the other ROI\n"
        "                        commands(-sphere, -cube)\n"
    );
    
    printf(
       "   -v[erbose]          : show debug information\n"
       "\n"
    );
    
    return 0;
}

int main(int argc, char * argv[])
{
    FSLIO *fslio;
	void *buffer;
	size_t buffsize;
	
	char *inputpath        = NULL;
	char *maskpath         = NULL;
	
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
    
    double sphereradius    = -1.0;
    FloatPoint3D center    = MakeFloatPoint3D(-9999.9, -9999.9, -9999.9);
    FloatPoint3D cubeDim   = MakeFloatPoint3D(-9999.9, -9999.9, -9999.9);
    BOOL resetVolume       = TRUE;
    BOOL verbose           = FALSE;
    long nSamples          = -1;
	
	int ac;
    
	if( argc < 2 ) return show_help();
    
    /* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if( ! strncmp(argv[ac], "-h", 2) ) {
			return show_help();
		} else if ( ! strcmp(argv[ac], "-input") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -input\n");
				return 1;
			}
			inputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-center") ) {
			if( ac+3 >= argc ) {
				fprintf(stderr, "** missing argument for -center, 3 coordinates must be supplied!\n");
				return 1;
			}
			ac++;
			float x = atoi(argv[ac]);
			ac++;
			float y = atoi(argv[ac]);
			ac++;
			float z = atoi(argv[ac]);
            center = MakeFloatPoint3D(x,y,z);
		} else if ( ! strcmp(argv[ac], "-cube") ) {
			if( ac+3 >= argc ) {
				fprintf(stderr, "** missing argument for -cube, 3 coordinates must be supplied!\n");
				return 1;
			}
			ac++;
			float x = atoi(argv[ac]);
			ac++;
			float y = atoi(argv[ac]);
			ac++;
			float z = atoi(argv[ac]);
            cubeDim = MakeFloatPoint3D(x,y,z);
		} else if ( ! strcmp(argv[ac], "-sphere") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -sphere\n");
				return 1;
			}
			sphereradius = atof(argv[ac]);
		} else if ( ! strcmp(argv[ac], "-randomsample") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -randomsample\n");
				return 1;
			}
			nSamples = atol(argv[ac]);
		} else if ( ! strncmp(argv[ac], "-m", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m\n");
				return 1;
			}
			maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-keepVolume") ) {
			resetVolume = FALSE;
		} else if ( ! strncmp(argv[ac], "-v", 2) ) {
			verbose = TRUE;
		} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
		}
	}
	
	if ( inputpath == NULL ) {
		fprintf(stderr, "No input volume specified(-input)!\n");
		return 1;
	}
	
	if ( maskpath == NULL ) {
		fprintf(stderr, "An output path for the mask be specified(-mask)!\n");
		return 1;
	}
    
    BOOL inputEqualsOutput = ! strcmp(inputpath, maskpath);
	
	if ( center.x < -9999.0 && nSamples < 0 ) {
		fprintf(stderr, "ROI center needs to be specified!(-center)!\n");
		return 1;
	}
    
    if ( sphereradius <= 0 && cubeDim.x < 0 && nSamples < 0 ) {
        fprintf(stderr, "ROI sphere radius or cube dimensions needs to be specified!(-sphere, -cube)!\n");
		return 1;
    }
	
    if ( verbose ) {
        fprintf(stdout, "Input file: %s\n", inputpath);
        fprintf(stdout, "Mask file: %s\n", maskpath);
        if ( sphereradius > 0 || cubeDim.x > -9999.9 )
        fprintf(stdout, "Center: %.2fmm %.2fmm %.2fmm\n", center.x, center.y, center.z);
        if ( sphereradius > 0) {
            fprintf(stdout, "Sphere radius: %.4fmm\n", sphereradius);
        }
        if ( cubeDim.x > -9999.9 ) {
            fprintf(stdout, "Cube: %.2fmm %.2fmm %.2fmm\n", cubeDim.x, cubeDim.y, cubeDim.z);
        }
    }
    
    fslio = FslOpen(inputpath, "rb");
    if (fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n",inputpath);
        return 1;
    }
    
	/* determine dimensions */
	FslGetDim(fslio, &xDim, &yDim, &zDim, &vDim);
    
    if ( verbose ) {
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", xDim, yDim, zDim, vDim);
    }
    
    if (fslio->niftiptr->scl_slope != 0) {
        slope = fslio->niftiptr->scl_slope;
        inter = fslio->niftiptr->scl_inter;
    }
	
	/* determine datatype and initalize buffer */
	dt = FslGetDataType(fslio, &pixtype);
	
    if ( verbose ) {
        fprintf(stdout, "Dt: %ld Pixtype: %d\n", dt, pixtype);
    }
    
    /* prepare mask file */
    FSLIO *fslioMask;
   	/* Init buffer */
    buffsize = (unsigned long)xDim * (unsigned long)yDim * (unsigned long)zDim * (unsigned long)(dt/8);
    buffer   = malloc(buffsize);
    
    /* Read in first volume */
    if (!FslReadVolumes(fslio, buffer, 1)) {
        free(buffer);
        fprintf(stderr, "\nError - reading data in %s\n", inputpath);
        FslClose(fslio);
        return 1;
    }
    
    double ***mask = FslGetVolumeAsScaledDouble(fslio, 0);
    
    
    /* Get sForm/qForm for MM<->Voxel space conversion */
    int sform_code, qform_code;
    mat44 sform44, qform44, stdmat44;
    
    sform_code = FslGetStdXform(fslio, &sform44);
    qform_code = FslGetRigidXform(fslio, &qform44);
    
    if (sform_code!=NIFTI_XFORM_UNKNOWN) {
        stdmat44 = sform44;
    } else if (qform_code!=NIFTI_XFORM_UNKNOWN) {
        stdmat44 = qform44;
    }
    
    if ( inputEqualsOutput ) {
        fprintf(stdout, "*Input equals output\n");
    } else {
        fprintf(stdout, "*Creating mask\n");
    }

    /* Preparing new mask file */
    FslClose(fslio);
    fslioMask = FslOpen(maskpath, "wb");
    FslCloneHeader(fslioMask, fslio);
    FslSetDim(fslioMask, xDim, yDim, zDim, 1);
    FslSetDimensionality(fslioMask, 4);
    FslSetDataType(fslioMask, pixtype);
    FslWriteHeader(fslioMask);
    free(fslio);
    
    if (fslioMask == NULL) {
        fprintf(stderr, "\nWarning, could not open %s for writing.\n", maskpath);
    }
    
    /* Iterate over all voxels */
    for (short z=0; z<zDim; z=z+1) {
        for (short y=0; y<yDim; y=y+1) {
            for (short x=0; x<xDim; x=x+1) {
                
                /* Convert current coordinate to MM space */
                float mmx = 0.0, mmy = 0.0, mmz = 0.0;
                FslGetMMCoord(stdmat44, x, y, z, &mmx, &mmy, &mmz);
                FloatPoint3D point = MakeFloatPoint3D(mmx, mmy, mmz);
                
                if ( sphereradius > 0 && rsVoxelInSphere(point, center, sphereradius) ) {
                    mask[z][y][x] = 1.0;
                } else if(cubeDim.x > -9999.9 && rsVoxelInCube(point, center, cubeDim)) {
                    mask[z][y][x] = 1.0;                    
                } else if ( resetVolume ) {
                    mask[z][y][x] = 0.0;
                }
            }
        }
    }
    
    /* If sampling is requested sample nSamples voxels from the mask */
    if (nSamples > 0) {
        if ( verbose ) fprintf(stdout, "Sampling %ld voxels\n", nSamples);
        
        // Create a new empty mask
        double ***oldMask = mask;
        mask = NULL;
        mask = d3matrix(zDim-1, yDim-1, xDim-1);
        
        for (short z=0; z<zDim; z=z+1) {
            for (short y=0; y<yDim; y=y+1) {
                for (short x=0; x<xDim; x=x+1) {
                    mask[z][y][x] = 0.0;
                }
            }
        }

        // Randomly take sample points and copy them from the original mask
        gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
        
        while (nSamples > 0) {
            unsigned short x = gsl_rng_uniform_int(r, xDim);
            unsigned short y = gsl_rng_uniform_int(r, yDim);
            unsigned short z = gsl_rng_uniform_int(r, zDim);
            
            if (oldMask[z][y][x] <= 0.5) {
                continue;
            }
            
            if (mask[z][y][x] > 0.5) {
                continue;
            }
            
            mask[z][y][x] = oldMask[z][y][x];
            nSamples = nSamples - 1L;
        }
        
        gsl_rng_free(r);
        free(oldMask);
    }
    
    convertScaledDoubleToBuffer(fslioMask->niftiptr->datatype, buffer, mask[0][0], slope, inter, xDim, yDim, zDim, FALSE);
    FslWriteVolumes(fslioMask, buffer, 1);
    free(buffer);
    free(mask);
    FslClose(fslioMask);
    free(fslioMask);
    
	return 0;
}