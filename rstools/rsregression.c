/*******************************************************************
 *
 * rstimecourse.c
 *
 * Extracts the time course from a 4D-Nifti by supplying either a voxel or a binary mask
 *
 * Usage: rstimecourse [-m <mask> [-a <algorithm>]] [-p <X> <Y> <Z>] -input <volume>
 *
 * 
 * André Hoffmann
 *******************************************************************/


#include <stdio.h>
#include <strings.h>

#include <nifti1.h>
#include <fslio.h>
#include "rsniftiutils.h"
#include "rsmathutils.h"

int show_help( void )
{
   printf(
      "rstimecourse: Given a 4D-Nifti, this tool extracts the time course\n"
      "              for a single voxel or the meaned average of a region\n"
      "              specified by a binary mask\n"
      "\n"
   );
    
   printf(
      "basic usage:  rstimecourse [-m <mask> [-a <algorithm>] [-savemask <mask>]] [-p <X> <Y> <Z>] -input <volume>\n"
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
      "              -v               : show debug information\n"
      "\n"
   );
    
   return 0;
}

void testRegression()
{
    int nSamples        = 10;
    double signal[]     = {1.9, 0.6, 0.5, 0.8, -0.4, -0.9, -0.7, -0.1, -1.7, -0.2};
    int nRegressors     = 2;
    /*double regressors[] = {-2.0, -1.0, -0.8, -0.3, 0.0, 0.5, 0.6, 0.7, 1.0, 1.2}; */
    double **regressors = d2matrix(nRegressors, nSamples);
    double r1[]         = {-2.0, -1.0, -0.8, -0.3, 0.0, 0.5, 0.6, 0.7, 1.0, 1.2};
    double r2[]         = {0.3, 0.1, 0.9, 0.1, 1.3, -1.3, -0.7,  0.4, -1.2, 0.8};
    regressors[0]       = r1;
    regressors[1]       = r2;
    double betas[]      = {0.0, 0.0};
    double residuals[]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double fitted[]     = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    rsLinearRegression(nSamples, signal, nRegressors, regressors, betas, residuals, fitted, TRUE);
    
    for ( int i = 0; i < nSamples; i++) {
        printf("Residual %d: %.4f\n", i+1, residuals[i]);
    }
    
    for ( int i = 0; i < nRegressors+1; i++) {
        printf("Beta %d: %.4f\n", i+1, betas[i]);
    }
    
    for ( int i = 0; i < nSamples; i++) {
        printf("Fitted %d: %.4f\n", i+1, fitted[i]);
    }
}

int main(int argc, char * argv[])
{
    FSLIO *fslio;
	void *buffer;
	size_t buffsize;
	
	char *inputpath = NULL;
	char *maskpath = NULL;
    char *savemaskpath = NULL;
    char *saveBetasPath = NULL;
    char *saveResidualsPath = NULL;
    char *saveFittedPath = NULL;
	
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
    
    BOOL verbose = FALSE;
	
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
		} else if ( ! strcmp(argv[ac], "-betas") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -betas\n");
				return 1;
			}
			saveBetasPath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-residuals") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -residuals\n");
				return 1;
			}
			saveResidualsPath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-fitted") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -fitted\n");
				return 1;
			}
			saveFittedPath = argv[ac];  /* no string copy, just pointer assignment */
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
		} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
		}
	}
	
	if ( inputpath == NULL ) {
		fprintf(stderr, "No input volume specified(-input)!\n");
		return 1;
	}
	
	if ( maskpath == NULL ) {
		fprintf(stderr, "A binary mask must be specified(-mask)!\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Input file: %s\n", inputpath);
        fprintf(stdout, "Mask file: %s\n", maskpath);
        fprintf(stdout, "Residuals file: %s\n", saveResidualsPath);
        fprintf(stdout, "Fitted file: %s\n", saveFittedPath);
        fprintf(stdout, "Betas file: %s\n", saveBetasPath);
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

    FSLIO *fslioResiduals;
    if ( saveResidualsPath != NULL ) {

        fslioResiduals = FslOpen(saveResidualsPath, "wb");
        
        if (fslioResiduals == NULL) {
            fprintf(stderr, "\nError, could not read header info for %s.\n",saveResidualsPath);
            return 1;
        }
        
        FslCloneHeader(fslioResiduals, fslio);
        FslSetDim(fslioResiduals, xDim, yDim, zDim, vDim);
        FslSetDimensionality(fslioResiduals, 4);
        FslSetDataType(fslioResiduals, pixtype);
        FslWriteHeader(fslioResiduals);
        
        // copy contents
        buffsize = xDim*yDim*zDim*vDim*dt/8;
        buffer = malloc(buffsize);
        FslReadVolumes(fslio, buffer, vDim);
        FslWriteVolumes(fslioResiduals, buffer, vDim);
        
        double ***v = d3matrix(1, 1, vDim);

        for (int t=0; t<vDim; t=t+1) {
            v[0][0][t] = vDim-t+100;
        }
        
        convertScaledDoubleToBuffer(fslioResiduals->niftiptr->datatype, buffer, v, slope, inter, vDim, 1, 1);
        free(v[0][0]);
        free(v[0]);
        free(v);
        
        WriteTimeSeries(fslioResiduals, buffer, 49, 49, 49, vDim);
        WriteTimeSeries(fslioResiduals, buffer, 50, 50, 50, vDim);
        WriteTimeSeries(fslioResiduals, buffer, 51, 51, 51, vDim);
        
        FslClose(fslioResiduals);
        free(fslioResiduals);
        free(buffer);
        return 0;
    }
 
    unsigned long nPoints = 0L;
    double ***mask = d3matrix(zDim, yDim, xDim);
    Point3D *maskPoints = ReadMask(maskpath, xDim, yDim, zDim, &nPoints, savemaskpath, fslio, mask);
    if ( maskPoints == NULL) {
        fprintf(stderr, "\nError: Mask invalid.\n");
        FslClose(fslio);
        return 1;
    }
    
    /* Create output volumes */
    
    
    /* Iterate over all voxels that are to be regressed */
    buffsize = vDim*dt/8;
    
    if ( verbose ) {
        fprintf(stdout, "Buffsize: %ld\n", buffsize);
    }
    
    buffer = malloc(buffsize);
    
    for ( unsigned long i=0; i<nPoints; i=i+1L) {
        Point3D p = maskPoints[i];
    
        /* read out timecourse */
/*        FslReadTimeSeries(fslio, buffer, p.x, p.y, p.z, vDim); */
        double sum = 0;
    }
    
    free(mask);
    free(maskPoints);
    free(buffer);
    FslClose(fslio);
    free(fslio);
    
	return 0;
}