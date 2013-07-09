//
//  rscorrelation.c
//  rstools
//
//  Created by André Hoffmann on 6/7/13.
//
//

#include "rscorrelation_common.h"

int show_help( void )
{
    rsCorrelationPrintHelp();
    return 0;
}

int main(int argc, char * argv[]) {
	
	if( argc < 2 ) return show_help();
	
    struct rsCorrelationParameters p = rsCorrelationLoadParams(argc, argv);

    // Prepare buffer
    size_t buffsize = p.vDim*p.dt/8;
    void *buffer = malloc(buffsize);

    /* Iterate over all voxels for which the correlation is to be computed */
    double signal[p.vDim];
    for (short z=0; z<p.zDim; z=z+1) {
        if (p.verbose) fprintf(stdout, "Computing slice Z%03hd/%03hd\n", z+1, p.zDim);
        for (short y=0; y<p.yDim; y=y+1) {
            for (short x=0; x<p.xDim; x=x+1) {
                                
                /* If it's not in the mask skip it to improve the performance */
                if (p.mask != NULL && p.mask[z][y][x] < 0.1) {
                    
                    /* set the value in the correlation file to NaN so that it is skipped in later processing steps */
                    p.correlation[z][y][x] = log(-1.0);

                    continue;
                }
                
                /* read out timecourse */
                FslReadTimeSeries(p.fslio, buffer, x, y, z, p.vDim);
                convertBufferToScaledDouble(signal, buffer, (long)p.vDim, p.slope, p.inter, p.fslio->niftiptr->datatype);
                
                /* compute correlation */
                p.correlation[z][y][x] = rsZCorrelation(signal, p.regressor, (size_t)p.nRegressorValues);
            }
        }
    }

    /* Write correlation file */
    rsCorrelationWriteCorrelationFile(&p);
    
    /* Close everything */
    rsCorrelationFree(&p);
    free(buffer);
    
	return 0;
}