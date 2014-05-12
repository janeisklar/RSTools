//
//  rscorrelation2.c
//  rstools
//
//  Created by André Hoffmann on 7/9/13.
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
    
    // Load file
    // Prepare buffer
    size_t buffsize = rsGetBufferSize(p.xDim, p.yDim, p.zDim, p.vDim, p.dt);
    void  *buffer   = malloc(buffsize);
    
    if (buffer == NULL) {
        fprintf(stdout, "Not enough free memory :-(\n");
        return 1;
    }
    
    FslReadVolumes(p.fslio, buffer, p.vDim);
    
    // Prepare buffer
    short x,y,z;
    double *timecourse;
    
    /* Iterate over all voxels for which the correlation is to be computed */
    #pragma omp parallel num_threads(p.threads) private(y,x,timecourse) shared(p,buffer)
    {
        #pragma omp for schedule(guided)
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
                    double *fullTimecourse = malloc(p.vDim*sizeof(double));
                    rsExtractTimecourseFromBuffer(p.fslio, fullTimecourse, buffer, p.slope, p.inter, MakePoint3D(x, y, z), p.xDim, p.yDim, p.zDim, p.vDim);
                    
                    /* add the defined delay to the regressor and adjust the timecourse */
                    double *regressor = &p.regressor[(short)fmax(0, -1*p.delay)];
                    timecourse = &fullTimecourse[(short)fmax(0, p.delay)];
                    const size_t regressorLength = p.nRegressorValues - fabs(p.delay);
                    
                    /* compute correlation */
					if ( p.monteCarloRepetitions > 0 ) {
						p.correlation[z][y][x] = rsMonteCarloZCorrelation(timecourse, regressor, regressorLength, p.monteCarloRepetitions, p.monteCarloSampleSize);
					} else if ( p.conversionMode == RSTOOLS_CORRELATION_CONVERSION_Z ) {
                        p.correlation[z][y][x] = rsZCorrelation(timecourse, regressor, regressorLength);
                    } else if ( p.conversionMode == RSTOOLS_CORRELATION_CONVERSION_NONE ) {
                        p.correlation[z][y][x] = rsCorrelation(timecourse, regressor, regressorLength);
                    } else if ( p.conversionMode == RSTOOLS_CORRELATION_CONVERSION_T ) {
                        p.correlation[z][y][x] = rsTCorrelation(timecourse, regressor, regressorLength);
                    }
                    
                    free(fullTimecourse);
                }
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