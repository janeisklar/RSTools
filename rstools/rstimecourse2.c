//
//  rstimecourse2.c
//  rstools
//
//  Created by André Hoffmann on 9/11/13.
//
//

#include <stdio.h>
#include <strings.h>
#include <omp.h>

#include <nifti1.h>
#include <fslio.h>
#include "rsniftiutils.h"

int show_help( void )
{
   printf(
      "rstimecourse2: Given a 4D-Nifti, this tool extracts the time course\n"
      "               for a single voxel or the meaned average of a region\n"
      "               specified by a binary mask\n"
      "\n"
   );
    
   printf(
      "basic usage:  rstimecourse2 [-m <mask> [-a <algorithm>] [-savemask <mask>]] [-p <X> <Y> <Z>] -input <volume>\n"
      "\n"
   );
    
   printf(
      "options:\n"
      "              -a <algorithm>   : the algorithm used to aggregate the data within\n"
      "                                 a ROI, e.g. mean\n"
   );

   printf(
      "              -help            : show this help\n"
   );
 
   printf(
      "              -input <volume>  : the volume from which the timecourse will be extracted\n"
   );
   
   printf(
      "              -mask <mask>     : a mask specifying the ROI\n"
   );
   
   printf(
      "              -p <X> <Y> <Z>   : speficies a voxel using nifti coordinates(0-based) from\n"
      "                                 which the timecourse is to be extracted\n"
   );
    
   printf(
      "              -savemask <mask> : optional path where the rescaled mask specified with -mask\n"
      "                                 will be saved. The saved file with have the same dimensions\n"
      "                                 as the input volume.\n"
   );
   
   printf(
      "              -v[erbose]       : show debug information\n"
      "\n"
   );
    
   return 0;
}

int main(int argc, char * argv[])
{
	FSLIO *fslio;
	void *buffer;
	size_t buffsize;
	
	char *inputpath = NULL;
	char *maskpath = NULL;
    char *savemaskpath = NULL;
	
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
    int threads = 1;
    
    BOOL verbose = FALSE;
	
	int ac;

	if( argc < 2 ) return show_help();   /* typing '-help' is sooo much work */

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
		} else if ( ! strncmp(argv[ac], "-m", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m\n");
				return 1;
			}
			maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-s", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -savemask\n");
				return 1;
			}
			savemaskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-v", 2) ) {
			verbose = TRUE;
		} else if ( ! strncmp(argv[ac], "-p", 2) ) {
			if( ac+3 >= argc ) {
				fprintf(stderr, "** missing argument for -p, 3 coordinates must be supplied!\n");
				return 1;
			}
			ac++;
			x = atoi(argv[ac]);
			ac++;
			y = atoi(argv[ac]);
			ac++;
			z = atoi(argv[ac]);
		} else if ( ! strcmp(argv[ac], "-threads") ) {
    		if( ++ac >= argc ) {
    			fprintf(stderr, "** missing argument for -threads\n");
    			return 1;
    		}
    		threads = atoi(argv[ac]);
    	} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
		}
	}
	
	if ( inputpath == NULL ) {
		fprintf(stderr, "No input volume specified!\n");
		return 1;
	}
	
	if ( maskpath == NULL && (x<0||y<0||z<0) ) {
		fprintf(stderr, "Either a binary mask or a voxel coordinate must be specified!\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Input file: %s\n", inputpath);
        fprintf(stdout, "Mask file: %s\n", maskpath);
        fprintf(stdout, "Voxel: %d %d %d\n", x, y, z);
    }
	
    fslio = FslOpen(inputpath, "rb");
    if (fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n",inputpath);
        return 1;
    }
    
	/* open nifti dataset header */
	/*fslio = FslReadHeader(inputpath);
	if (fslio == NULL) {
		fprintf(stderr, "\nError, could not read header info for %s.\n",inputpath);
		return 1;
	}*/

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
	
    if ( x >= 0 ) {
        /*** Extract timecourse for a single point ***/
        
        /* check inputs */
        if ( (x<0) || (x>=fslio->niftiptr->nx) ) {
            fprintf(stderr, "\nError: x index (%d) out of range [0..%d]\n",x,fslio->niftiptr->nx-1);
            FslClose(fslio);
            return 1;
        }
        if ( (y<0) || (y>=fslio->niftiptr->ny) ) {
            fprintf(stderr, "\nError: y index (%d) out of range [0..%d]\n",y,fslio->niftiptr->ny-1);
            FslClose(fslio);
            return 1;
        }
        if ( (z<0) || (z>=fslio->niftiptr->nz) ) {
            fprintf(stderr, "\nError: z index (%d) out of range [0..%d]\n",z,fslio->niftiptr->nz-1);
            FslClose(fslio);
            return 1;
        }
        
        buffsize = vDim*dt/8;
        
        if ( verbose ) {
            fprintf(stdout, "Buffsize: %ld\n", buffsize);
        }
        
        buffer = malloc(buffsize);
        
        /* read out timecourse */
        FslReadTimeSeries(fslio, buffer, x, y, z, vDim);
        
        for (t=0; t<vDim; t++) {
            void *newV = (char*)buffer + (t*dt/8);
            
            double v = 0.0;
            int ret = convertBufferToScaledDouble(&v, newV, 1, slope, inter, fslio->niftiptr->datatype);
         
            if ( verbose ) {
                fprintf(stdout, "%d: %.2f(%d)\n", t, v, ret);
            } else {
                fprintf(stdout, "%.2f\n", v);
            }
        }
        
        /* clear buffer */
        free(buffer);
    } else {
        /*** Extract timecourse for a mask ***/
        unsigned long nPoints = 0L;
        Point3D *maskPoints = ReadMask(maskpath, xDim, yDim, zDim, &nPoints, savemaskpath, fslio, NULL);
        if ( maskPoints == NULL) {
            fprintf(stderr, "\nError: Mask invalid.\n");
            FslClose(fslio);
            return 1;
        }
        
        double timecourse[vDim];
        
        buffsize = (size_t)xDim*(size_t)yDim*(size_t)zDim*(size_t)vDim*(size_t)dt/(size_t)8;
        buffer = malloc(buffsize);
        FslReadVolumes(fslio, buffer, vDim);
        
        /* Read in volume */
        double ****volume = d4matrix(vDim-1, zDim-1, yDim-1, xDim-1);
        convertBufferToScaledDouble(volume[0][0][0], buffer, (long)xDim*(long)yDim*(long)zDim*(long)vDim, slope, inter, fslio->niftiptr->datatype);
        free(buffer);
        
        int t;

        /* Iterate over all timepoints */
        #pragma omp parallel num_threads(threads) private(t) shared(timecourse,nPoints)
        {
            #pragma omp for schedule(guided)
            for ( t = 0; t<vDim; t=t+1 ) {
            
                double sum = 0;
                
                /* Iterate over all points in the mask */
                for ( unsigned long i=0; i<nPoints; i=i+1) {
                    Point3D p = maskPoints[i];
                    sum = sum + volume[t][p.z][p.y][p.x];
                }
                
                /* Create average */
                timecourse[t] = sum / nPoints;
            }
        }
        
        free(volume[0][0][0]);
        free(volume[0][0]);
        free(volume[0]);
        free(volume);
        free(maskPoints);
        
        for ( t = 0; t<vDim; t=t+1 ) {
            fprintf(stdout, "%.10f\n", timecourse[t]);
        }
    }
    
    FslClose(fslio);
    free(fslio);
    
	return 0;
}
