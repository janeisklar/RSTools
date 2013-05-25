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

int main(int argc, char * argv[])
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
    
	return 0;
}
